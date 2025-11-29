#!/usr/bin/env python

from __future__ import print_function
from optparse import OptionParser
import json
import os
import pdb
import sys

import h5py
import numpy as np
import pandas as pd

import tensorflow as tf

import dataset
import seqnn
import trainer
import metrics
import math
import gc

"""
cello_train_head.py

Train a new dnn-head on top of an existing Basenji model.
"""


def make_pairs(boundaries, n_samples):
    pairs = []
    end_pairs = []
    last_b = -1
    for b in boundaries:
        for i in range(last_b + 1, b):
            if i + 1 < b:
                pairs.append((i, i + 1))
            if i + 1 == b:
                end_pairs.append((i, i + 1))
        last_b = b
    t_indices, t1_indices = zip(*pairs)
    v_indices, v1_indices = zip(*end_pairs)
    return (
        tf.convert_to_tensor(list(t_indices), dtype=tf.int32),
        tf.convert_to_tensor(list(t1_indices), dtype=tf.int32),
        tf.convert_to_tensor(list(v_indices), dtype=tf.int32),
        tf.convert_to_tensor(list(v1_indices), dtype=tf.int32),
    )


def rank_regression_loss(y_true, y_pred, temperature=1e-2):
    #
    point_loss = tf.reduce_mean(tf.square(y_true - y_pred))
    y_true = tf.math.log1p(y_true)
    y_pred = tf.math.log1p(y_pred)
    return tf.reduce_mean(tf.square(y_true - y_pred))


class MemoryClearCallback(tf.keras.callbacks.Callback):
    def on_epoch_end(self, epoch, logs=None):
        gc.collect()
        tf.keras.backend.clear_session()


@tf.function
def map_fn(
    idx,
    seq,
    fea,
    tar,
    mode,
    target_batch,
    train_indices=None,
    train1_indices=None,
    val_indices=None,
    val1_indicies=None,
):
    # Repeat sequences to match target batch size
    seq = tf.tile(seq, [target_batch, 1, 1])
    fea = tf.transpose(fea, perm=[1, 2, 0])
    fea = tf.squeeze(fea, axis=-1)
    tar = tf.transpose(tar, perm=[2, 1, 0])
    tar = tf.divide(tar, 128.0)
    if mode == "train":
        seq = tf.gather(seq, train_indices, axis=0)
        fea = tf.gather(fea, train_indices, axis=0)
        tar0 = tf.gather(tar, train_indices, axis=0)  # shape: (batch, 896, 1)
        tar1 = tf.gather(tar, train1_indices, axis=0)  # shape: (batch, 896, 1)
        tar = tf.stack([tar0, tar1], axis=1)
        return (seq, fea), tar

    if mode == "valid":
        seq = tf.gather(seq, val_indices, axis=0)
        fea = tf.gather(fea, val_indices, axis=0)
        tar0 = tf.gather(tar, val_indices, axis=0)  # shape: (batch, 896, 1)
        tar1 = tf.gather(tar, val1_indices, axis=0)  # shape: (batch, 896, 1)
        tar = tf.stack([tar0, tar1], axis=1)
        return (seq, fea), tar


@tf.function
def map_batch(idx, seq, fea, tar, mode, target_batch, train_indices=None, val_indices=None, shuffle=True):
    batch_size = tf.shape(seq)[0]

    def call_map_fn(i):
        return map_fn(i, seq[i : i + 1], fea[i : i + 1], tar[i : i + 1], mode, target_batch, train_indices, val_indices)

    (seqs, feas), tars = tf.map_fn(
        call_map_fn, tf.range(batch_size), fn_output_signature=((tf.float16, tf.float16), tf.float16)
    )

    # Concatenate outputs along batch dimension (0)
    seqs_flat = tf.reshape(seqs, [-1, 896, 1536])
    feas_flat = tf.reshape(feas, [-1, feature_len])  # Use 'feas' from map_fn output
    tars_flat = tf.reshape(tars, [-1, 896, 1])  # Use 'tars' from map_fn output

    return (seqs_flat, feas_flat), tars_flat


def feature_bytes(values):
    values = values.flatten().tobytes()
    return tf.train.Feature(bytes_list=tf.train.BytesList(value=[values]))


################################################################################
# main
################################################################################
def max_swish(x, clip=1.0):
    return tf.minimum(x * tf.nn.sigmoid(x), clip)


def stochastic_mask(x, rate=0.3):
    mask = tf.cast(tf.random.uniform(tf.shape(x)) > rate, tf.float32)
    return x * mask


def make_model(target_size=896, filters=1536, feature_len=1810, model_png="model.png"):
    from tensorflow.keras import layers, Model, regularizers

    # dense_pathway = feature_counts
    # dense_pathway = tf.keras.layers.GaussianNoise(stddev=0.2)(dense_pathway)
    #

    # Inputs
    sequence = tf.keras.Input(shape=(target_size, filters), name="sequence")
    feature_counts = tf.keras.Input(shape=(feature_len,), name="feature_counts")

    l2 = tf.keras.regularizers.l1_l2(l1=0, l2=0)
    dense_pathway = tf.keras.layers.Lambda(
        lambda x: tf.keras.backend.in_train_phase(x * tf.random.normal(tf.shape(x), mean=1.0, stddev=0.02), x)
    )(feature_counts)
    dense_pathway = tf.keras.layers.Lambda(lambda x: tf.math.log1p(x))(dense_pathway)
    # dense_pathway = tf.keras.layers.Lambda(lambda x: tf.keras.backend.in_train_phase(stochastic_mask(x, rate=0.50), x))(dense_pathway)
    denses = [128, 128]
    # dense_pathway = layers.Dropout(0.20)(dense_pathway)
    for i, w in enumerate(denses):
        # dense_pathway = layers.LayerNormalization()(dense_pathway)

        if i < 1:
            dense_pathway = layers.Dense(w, kernel_regularizer=l2, activation="elu")(dense_pathway)
            # dense_pathway = layers.Dropout(0.25)(dense_pathway)
        else:
            dense_pathway = layers.Dense(w, kernel_regularizer=l2, activation="elu")(dense_pathway)
            pass

    # gate = tf.keras.layers.Dense(filters, activation="linear")(dense_pathway)

    gate = layers.Dense(
        filters,
        activation="linear",
        kernel_initializer=tf.keras.initializers.RandomNormal(mean=0.0, stddev=0.01),
        bias_initializer=tf.keras.initializers.RandomNormal(0.0, 0.01),  # Start near 1.0
    )(dense_pathway)
    gate = tf.keras.layers.RepeatVector(target_size)(gate)
    current = layers.Multiply()([gate, sequence])

    # current = layers.Add()([sequence, current])
    # current = layers.Conv1D(filters=128, kernel_size=3, padding='same', activation='swish')(current)

    # onoff = tf.keras.layers.Dense(filters, activation="gelu", kernel_initializer='he_normal')(dense_pathway)
    # onoff = tf.keras.layers.RepeatVector(target_size)(onoff)
    # current = tf.keras.layers.Concatenate(axis=-1)([current, onoff])

    # gate2 = tf.keras.layers.Dense(filters, activation="elu", kernel_initializer='he_normal')(dense_pathway)
    # gate2 = tf.keras.layers.RepeatVector(target_size)(gate2)
    # current = tf.keras.layers.Concatenate(axis=-1)([current, gate2])

    # current = tf.keras.layers.TimeDistributed(tf.keras.layers.Dense(512, activation='elu'))(current)
    # current = tf.keras.layers.TimeDistributed(tf.keras.layers.Dense(64, activation='relu'))(current)
    # current = tf.keras.layers.TimeDistributed(tf.keras.layers.Dense(1, activation='softplus'))(current)
    current = tf.reduce_sum(current, axis=-1, keepdims=True)
    current = tf.keras.layers.Activation("softplus")(current)
    # Build and compile model
    model = Model(inputs=[sequence, feature_counts], outputs=current)
    model.compile()

    # Plot model architecture
    tf.keras.utils.plot_model(model, to_file=model_png, show_shapes=True)
    print(model.summary())

    return model


def main():
    usage = "usage: %prog [options] <params_file> <model_file> <data_dir>"
    parser = OptionParser(usage)
    parser.add_option(
        "-o", dest="out_dir", default="head_out", help="Output directory for new head training [Default: %default]"
    )
    parser.add_option("--shifts", dest="shifts", default="0", help="Ensemble prediction shifts [Default: %default]")
    (options, args) = parser.parse_args()

    if len(args) != 3:
        parser.error("Must provide parameters, model, and test data directory")
    else:
        params_file = args[0]
        model_file = args[1]
        data_dir = args[2]

    os.makedirs(options.out_dir, exist_ok=True)
    head_dir = f"{data_dir}_head"
    os.makedirs(head_dir, exist_ok=True)
    # parse shifts to integers
    options.shifts = [int(shift) for shift in options.shifts.split(",")]

    #######################################################
    # prepare model

    # read model parameters
    with open(params_file) as params_open:
        params = json.load(params_open)
    params_model = params["model"]
    params_train = params["train"]

    #######################################################
    # predict
    target_batch = 63  # 113
    if os.path.exists(f"{head_dir}/tfrecords"):
        print("TF Records Exist!")
    else:
        with open(f"{data_dir}/statistics.json") as stats_open:
            stats = json.load(stats_open)
        stats[
            "num_targets"
        ] = target_batch  # red update data_dir stats to have number of targets equal to tissue samples
        with open(f"{data_dir}/statistics.json", "w") as stats_open:
            json.dump(stats, stats_open)
        stats["seq_depth"] = 1536
        stats["seq_1hot"] = False
        stats["seq_length"] = 896
        with open(f"{head_dir}/statistics.json", "w") as stats_open:
            json.dump(stats, stats_open)
        os.makedirs(f"{head_dir}/tfrecords", exist_ok=True)
        seqnn_model = seqnn.SeqNN(params_model)
        seqnn_model.restore(model_file, trunk=True)
        seqnn_model.build_embed(-1)
        seqnn_model.build_ensemble(False, options.shifts)

        for split_label in ["train", "valid", "test"]:
            split_data = dataset.SeqDataset(
                data_dir, split_label=split_label, batch_size=params_train["batch_size"], mode="eval"
            )
            print(split_data.batches_per_epoch())
            split_preds = []
            split_feats = []
            split_targets = []

            batch_threshold = 1024

            split_iter = iter(split_data.dataset)
            keep_count = 0
            tfr_num = 0
            for si in range(split_data.batches_per_epoch()):
                keep_count += 1
                x, x1, y = next(split_iter)

                preds = seqnn_model.predict(x, verbose=0)
                split_preds.append(preds)
                split_feats.append(x1)
                split_targets.append(y)

                if len(split_feats) == batch_threshold or si == split_data.batches_per_epoch() - 1:
                    split_preds = np.vstack(split_preds)
                    split_feats = np.vstack(split_feats)
                    split_targets = np.vstack(split_targets)
                    tf_opts = tf.io.TFRecordOptions(compression_type="ZLIB")
                    print(f"Writing {split_label}-{tfr_num}.tfr")
                    with tf.io.TFRecordWriter(f"{head_dir}/tfrecords/{split_label}-{tfr_num}.tfr", tf_opts) as writer:
                        for p, f, t in zip(split_preds, split_feats, split_targets):
                            # print(p.shape, f.shape, t.shape)
                            p = p.astype(np.float16)
                            f = f.astype(np.float16)
                            t = t.astype(np.float16)
                            features_dict = {
                                "sequence": feature_bytes(p),
                                "feature_counts": feature_bytes(f),
                                "target": feature_bytes(t),
                            }
                            example = tf.train.Example(features=tf.train.Features(feature=features_dict))
                            writer.write(example.SerializeToString())
                    tfr_num += 1

                    split_preds = []
                    split_feats = []
                    split_targets = []

    #######################################################
    # Training

    global feature_len
    feature_len = 1810  # 1575
    train_data = [
        dataset.SeqDataset(
            head_dir,
            split_label="train",
            batch_size=1,
            shuffle_buffer=params_train.get("shuffle_buffer", 128),
            mode="train",
        )
    ]
    eval_data = [dataset.SeqDataset(head_dir, split_label="valid", batch_size=1, mode="valid")]

    model = make_model(model_png=f"{options.out_dir}/model.png")

    num_targets = model.output_shape[-1]
    model_metrics = [
        metrics.PearsonR(num_targets),
        metrics.R2(num_targets),
        metrics.SeqAUC(curve="PR"),
        metrics.SeqAUC(curve="ROC"),
    ]

    # 0.002
    """
    initial_learning_rate = 0.002
    optimizer = tf.keras.optimizers.AdamW(
        learning_rate=initial_learning_rate,
        weight_decay=0.002,
        beta_1=0.9,
        beta_2=0.999,
        clipnorm=10,
    )
    """
    initial_learning_rate = 0.002
    optimizer = tf.keras.optimizers.SGD(
        learning_rate=initial_learning_rate,
        momentum=0.9,  # optional but recommended
        nesterov=True,  # optional: adds a slight boost in many cases
        clipnorm=0.5,  # gradient clipping
    )

    checkpoint_callback = tf.keras.callbacks.ModelCheckpoint(
        f"{options.out_dir}/best_model.h5", save_best_only=True, mode="max", monitor="val_pearsonr", verbose=1
    )
    memory_clear_callback = MemoryClearCallback()

    # 38, 40, 41
    # train_exempt = [6, 12, 19, 26, 32, 39, 44, 50, 53, 56, 59, 62]  # [5, 12, 18, 25, 32, 38, 44, 49, 52, 55, 58, 61]
    # train_indices = list(set(range(target_batch)) - set(train_exempt))
    # train_indices = tf.convert_to_tensor(train_indices, dtype=tf.int32)
    # val_indices = list(set([6, 12, 19, 26, 32, 39, 44, 50, 53, 56, 59, 62]))
    # val_indices = tf.convert_to_tensor(val_indices, dtype=tf.int32)

    boundaries = [6, 12, 19, 26, 32, 39, 44, 50, 53, 56, 59, 62]  # example boundaries
    n_samples = 63
    t_indices, t1_indices, v_indices, v1_indices = make_pairs(boundaries, n_samples)
    train_first = eval_data[0].dataset
    iterator = iter(train_first)
    data = next(iterator)
    # seq, fea, tar

    print(tf.shape(data[2]))
    # sys.exit()
    x, y = map_fn(0, *data, "train", target_batch, train_indices, val_indices)
    np.set_printoptions(formatter={"float_kind": "{:.6f}".format})
    print(x[1].numpy()[:, 1182])
    print(tf.shape(x[0]), tf.shape(x[1]), tf.shape(y))
    sys.exit()
    lr_scheduler = tf.keras.callbacks.ReduceLROnPlateau(
        monitor="val_pearsonr", factor=0.5, patience=3, verbose=1, min_lr=1e-6, mode="max"
    )

    # Add index to each example in the dataset
    train_ds = (
        train_data[0]
        .dataset.enumerate()
        .map(
            lambda idx, data: map_batch(idx, *data, "train", target_batch, train_indices, val_indices),
            num_parallel_calls=4,
        )
    )
    valid_ds = (
        eval_data[0]
        .dataset.enumerate()
        .map(
            lambda idx, data: map_batch(idx, *data, "valid", target_batch, train_indices, val_indices),
            num_parallel_calls=4,
        )
        .repeat()
    )
    print(eval_data[0].batches_per_epoch())
    """
    model.compile(loss=pointwise, optimizer=optimizer, metrics=model_metrics)
    model.fit(
        train_ds,
        epochs=1,
        steps_per_epoch=1500,#train_data[0].batches_per_epoch(),
        validation_data=valid_ds,
        validation_steps=eval_data[0].batches_per_epoch(),
        callbacks=[checkpoint_callback, memory_clear_callback, lr_scheduler],
        # verbose=2,
    )
    """
    model.compile(loss=rank_regression_loss, optimizer=optimizer, metrics=model_metrics)
    for epoch_big in range(50):
        print(epoch_big)
        model.fit(
            train_ds,
            epochs=1,
            steps_per_epoch=5000,  # train_data[0].batches_per_epoch(),
            validation_data=valid_ds,
            validation_steps=eval_data[0].batches_per_epoch(),
            callbacks=[checkpoint_callback, lr_scheduler],
            # verbose=2,
        )
        np.set_printoptions(suppress=True, precision=2)
        for ds in (train_ds, valid_ds):
            x_batch, y_batch = next(iter(ds))
            preds = model.predict(x_batch)
            mean_true = tf.math.reduce_mean(tf.squeeze(y_batch, axis=-1), axis=1)
            mean_pred = tf.math.reduce_mean(tf.squeeze(preds, axis=-1), axis=1)
            z_true = (mean_true.numpy() - mean_true.numpy().mean()) / mean_true.numpy().std()
            z_pred = (mean_pred.numpy() - mean_pred.numpy().mean()) / mean_pred.numpy().std()
            print(mean_pred.numpy())
            print(mean_true.numpy())
            print(z_true - z_pred)
        # print(z_pred)

    #######################################################
    # Evaluation


################################################################################
# __main__
################################################################################
# if __name__ == "__main__":
#    main()
