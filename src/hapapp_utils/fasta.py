"""Reading FASTA files and fetching subsequences, with no index file written next to them."""

from __future__ import annotations

from collections import defaultdict


def iter_fasta(path):
    """Yield (name, sequence) for each record. The name is the first word of the header."""
    name = None
    chunks = []
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def read_lengths(path):
    """Sequence name -> length, without holding the sequences in memory."""
    lengths = {}
    name = None
    with open(path) as handle:
        for line in handle:
            if line.startswith(">"):
                name = line[1:].split()[0]
                lengths[name] = 0
            elif name is not None:
                lengths[name] += len(line.strip())
    return lengths


def count_records(path):
    with open(path) as handle:
        return sum(1 for line in handle if line.startswith(">"))


def write_fasta(records, path):
    with open(path, "w") as out:
        for name, sequence in records:
            out.write(f">{name}\n{sequence}\n")


def fetch_subsequences(fasta_path, keys_path, out_path, rename=None):
    """Write the subsequences listed in a keys file, in the order listed.

    Each key line is: new name, from, to, source name. Coordinates are 1-based and inclusive, as
    for esl-sfetch -Cf. A range running past an end of its source is cut short, and the caller can
    judge the result by its length. rename maps source names in the keys to names in the FASTA.
    """
    rename = rename or {}
    keys = []
    wanted = defaultdict(list)
    with open(keys_path) as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            new_name, start, end, source = line.split()[:4]
            source = rename.get(source, source)
            wanted[source].append((len(keys), int(start), int(end)))
            keys.append(new_name)

    sequences = [None] * len(keys)
    for name, sequence in iter_fasta(fasta_path):
        for index, start, end in wanted.pop(name, ()):
            sequences[index] = sequence[max(start, 1) - 1:end]
        if not wanted:
            break
    if wanted:
        missing = sorted(wanted)
        raise ValueError(f"{len(missing)} sequence name(s) in the keys are not in {fasta_path}: {missing[:5]}")
    write_fasta(zip(keys, sequences), out_path)
    return len(keys)
