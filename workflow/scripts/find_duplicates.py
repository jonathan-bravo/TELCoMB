#!/usr/bin/env python

from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-p', help='BLAT alignment file (.psl)', required=True)
    parser.add_argument('-o', required=True)
    parser.add_argument('-s', help='Similarity threshold', type=float, required=True)
    return parser.parse_args()


def read_psl(psl):
    # return [q.strip().split('\t') for q in open(psl)]
    return (q.strip().split('\t') for q in open(psl))


def skip_header(qresults):
    try: [next(qresults) for _ in range(5)]
    except StopIteration: pass


def high_match(q, threshold): 
    return (int(q[0]) >= threshold * int(q[14]) # amount covered is 90% of target?
    and int(q[0]) >= threshold * int(q[10]) # amount covered is 90% of query?
    and q[9] != q[13]) # doesn't equal self


def find_dupes(qresults, threshold):
    dupes = {}
    [dupes.setdefault(q[9], []).append(q[13]) for q in qresults if high_match(q, threshold)]
    return dupes


def try_remove(dupes, x):
    try: dupes.pop(x)
    except KeyError: pass


def look_through_dupes(dupes, k):
    try: [try_remove(dupes, x) for x in dupes[k]]
    except KeyError: pass


def remove_doubles(dupes):
    results = []
    [look_through_dupes(dupes, k) for k in list(dupes)]
    [results.extend(dupes[s]) for s in dupes]
    return set(results)


def write_out(outfile, dup_set):
    with open(outfile, 'w') as o:
        [o.write(f'{dup}\n') for dup in dup_set]


def main():
    args = parse_args()
    qresults = read_psl(args.p)
    skip_header(qresults)
    dupes = find_dupes(qresults, args.s)
    cleaned_dupes = remove_doubles(dupes)
    write_out(args.o, cleaned_dupes)


if __name__ == "__main__":
    main()