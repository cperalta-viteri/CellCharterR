#!/usr/bin/env python3
"""Export a pinned PyCave GaussianMixture state dict without importing torch.

The official CellCharter AutoK archive stores a small PyTorch ZIP/pickle state
dict. This exporter decodes only the tensor/storage opcodes used by that state
dict and writes portable JSON for one requested K. It is a fixture-generation
utility, not a CellCharterR runtime dependency.
"""

import argparse
import io
import json
import pickle
import struct
import zipfile


class _FloatStorage:
    format = "f"
    width = 4


class _DoubleStorage:
    format = "d"
    width = 8


class _Storage:
    def __init__(self, values):
        self.values = values


class _Tensor:
    def __init__(self, storage, offset, size, stride):
        self.storage = storage
        self.offset = int(offset)
        self.size = tuple(int(value) for value in size)
        self.stride = tuple(int(value) for value in stride)

    def to_list(self):
        if not self.size:
            return self.storage.values[self.offset]

        def value_at(index):
            position = self.offset + sum(i * step for i, step in zip(index, self.stride))
            return self.storage.values[position]

        def build(prefix, axis):
            if axis == len(self.size):
                return value_at(prefix)
            return [build(prefix + (i,), axis + 1) for i in range(self.size[axis])]

        return build(tuple(), 0)


def _rebuild_tensor_v2(storage, offset, size, stride, _requires_grad, _hooks):
    return _Tensor(storage, offset, size, stride)


class _TorchUnpickler(pickle.Unpickler):
    def __init__(self, stream, archive, prefix):
        super().__init__(stream)
        self.archive = archive
        self.prefix = prefix

    def find_class(self, module, name):
        if module == "torch._utils" and name == "_rebuild_tensor_v2":
            return _rebuild_tensor_v2
        if module == "torch" and name == "FloatStorage":
            return _FloatStorage
        if module == "torch" and name == "DoubleStorage":
            return _DoubleStorage
        return super().find_class(module, name)

    def persistent_load(self, identifier):
        if not isinstance(identifier, tuple) or identifier[0] != "storage":
            raise pickle.UnpicklingError(f"Unsupported persistent identifier: {identifier!r}")
        storage_type, key, size = identifier[1], str(identifier[2]), int(identifier[4])
        raw = self.archive.read(f"{self.prefix}data/{key}")
        expected = size * storage_type.width
        if len(raw) != expected:
            raise ValueError(f"Storage {key} has {len(raw)} bytes; expected {expected}.")
        values = struct.unpack(f"<{size}{storage_type.format}", raw)
        return _Storage(values)


def export_state(outer_zip, k):
    suffix = f"/GaussianMixture_k{k}/model/parameters.pt"
    with zipfile.ZipFile(outer_zip) as archive:
        matches = [name for name in archive.namelist() if name.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one {suffix} entry, found {len(matches)}.")
        parameter_bytes = archive.read(matches[0])

    with zipfile.ZipFile(io.BytesIO(parameter_bytes)) as archive:
        data_pickle = next(name for name in archive.namelist() if name.endswith("data.pkl"))
        prefix = data_pickle[: -len("data.pkl")]
        state = _TorchUnpickler(io.BytesIO(archive.read(data_pickle)), archive, prefix).load()

    return {
        "format": "PyCave GaussianMixture state_dict exported from parameters.pt",
        "k": k,
        "component_probs": state["component_probs"].to_list(),
        "means": state["means"].to_list(),
        "precisions_cholesky": state["precisions_cholesky"].to_list(),
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("archive", help="Official CellCharter AutoK ZIP")
    parser.add_argument("output", help="Destination JSON")
    parser.add_argument("--k", type=int, default=11)
    args = parser.parse_args()
    state = export_state(args.archive, args.k)
    with open(args.output, "w", encoding="utf-8") as stream:
        json.dump(state, stream, separators=(",", ":"))
        stream.write("\n")


if __name__ == "__main__":
    main()
