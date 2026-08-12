#!/usr/bin/env python3

import argparse

import numpy as np
from sklearn.metrics import pair_confusion_matrix


def load_labels(label_file):
    labels = np.loadtxt(label_file, dtype=np.int64, ndmin=2)
    if labels.shape[0] != 2:
        raise ValueError(
            f"{label_file} must contain exactly two rows: "
            "predicted labels followed by ground-truth labels"
        )
    if labels.shape[1] == 0:
        raise ValueError(f"{label_file} does not contain any labels")
    return labels


def compute_metrics(label_file):
    labels = load_labels(label_file)
    y_pred = labels[0, :]
    y_true = labels[1, :]

    conf_matrix = pair_confusion_matrix(y_true, y_pred)
    a, b = conf_matrix[0]
    c, d = conf_matrix[1]

    precision = d / (b + d) if (b + d) > 0 else 0
    recall = d / (c + d) if (c + d) > 0 else 0
    f1 = 2 * precision * recall / (precision + recall) if (precision + recall) > 0 else 0

    print(f"{label_file}")
    print(f"Precision: {precision:.4f}")
    print(f"Recall:    {recall:.4f}")
    print(f"F1-score:  {f1:.4f}")


def main():
    parser = argparse.ArgumentParser(
        description="Compute pairwise precision, recall, and F1 from a two-row label file."
    )
    parser.add_argument(
        "label_file",
        help="two-row file containing predicted labels and ground-truth labels",
    )
    args = parser.parse_args()
    compute_metrics(args.label_file)


if __name__ == "__main__":
    main()
