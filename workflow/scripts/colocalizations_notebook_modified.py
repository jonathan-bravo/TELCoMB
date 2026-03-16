#!/usr/bin/env python

import seaborn as sns
import json
from matplotlib import pyplot as plt
from matplotlib.patches import Rectangle
from argparse import ArgumentParser


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('--read_lengths')
    parser.add_argument('--colocalizations')
    parser.add_argument('--output_plot')
    return parser.parse_args()
    

def add_genes_to_plot(gene, position, color, label, axis):
    length = int(position.split(':')[1]) - int(position.split(':')[0])
    start = int(position.split(':')[0])
    element_r = Rectangle(
        (start, curr_y),
        length,
        height,
        facecolor=color,
        alpha=0.8,
        edgecolor='black',
        label=label
    )
    axis.add_patch(element_r)

    rx, ry = element_r.get_xy()
    cx = rx + element_r.get_width() / 2.0
    cy = ry + element_r.get_height() / 2.0
    axis.annotate(
        gene,
        (cx, cy - (height / 2 + 1)),
        color='black',
        fontsize=8,
        ha='center',
        va='center'
    )


def get_read_lengths(read_lengths):
    with open(read_lengths, 'r') as json_file:
        return json.load(json_file)
    

def get_colocalizations(colocalizations):
    colocalizations_dict = {}
    
    with open(colocalizations, 'rt') as c:
        next(c)
        for row in c:
            items = row.strip().split(',')
            read_id = items[0]
            args = items[1].split(';')
            args_pos = items[2].split(';')
            mges = items[3].split(';')
            mges_pos = items[4].split(';')
            for i, arg in enumerate(args):
                for j, mge in enumerate(mges):
                    a = arg.split("|")[4]
                    m = ':'.join(mge.split(':')[1:])
                    try:
                        colocalizations_dict[f'{a},{m}']['count'] += 1
                    except:
                        colocalizations_dict[f'{a},{m}'] = {
                            'read_id': read_id,
                            'arg_pos': args_pos[i],
                            'mge_pos': mges_pos[j],
                            'count': 1,
                        }
    return colocalizations_dict


def main():
    global curr_y
    global height

    ARG_COLOR  = 'red'
    MGE_COLOR  = 'green'

    args = parse_args()

    output_plot = args.output_plot

    read_lengths = get_read_lengths(args.read_lengths)
    colocalizations = get_colocalizations(args.colocalizations)

    print(len(colocalizations))    

    if len(colocalizations) > 0:
        used_reads_lengths = list()
        padding = 3
        height = 2
        curr_y = 1 + padding

        # Calculate the maximum count width for proper spacing
        max_count = max(item['count'] for item in colocalizations.values())
        count_width = len(str(max_count)) * 8  # Approximate width per digit
        left_margin = count_width + 30  # Add some padding

        plt.figure(figsize=(16, 3*(len(colocalizations)/7)))
        currentAxis = plt.gca()
        for item in colocalizations:
            arg, mge = item.split(',')
            read_length = read_lengths[colocalizations[item]['read_id']]
            used_reads_lengths.append(read_length)
            read_r = Rectangle(
                (0, curr_y+0.75),
                read_length,
                0.5,
                facecolor='black',
                alpha=0.8,
                edgecolor='black',
                label='Read',
            )
            currentAxis.add_patch(read_r)
            currentAxis.text(
                -left_margin/100,
                curr_y + height / 2,
                f"{colocalizations[item]['count']} X", 
                verticalalignment='center',
                horizontalalignment='right',
                fontsize=8,
                color='black'
            )
            add_genes_to_plot(
                arg,
                colocalizations[item]['arg_pos'],
                ARG_COLOR,
                'ARG',
                currentAxis
            )
            add_genes_to_plot(
                mge,
                colocalizations[item]['mge_pos'],
                MGE_COLOR,
                'MGE',
                currentAxis
            )

            curr_y += height + padding

        sns.despine(top=True, right=True, left=True, bottom=False)
        plt.tick_params(
            axis='y',
            which='both',
            left=False,
            right=False,
            labelleft=False
        )
        plt.xlim([-left_margin/100, max(used_reads_lengths)])
        plt.ylim([0, curr_y])
        # plt.grid(axis='x', color='0.95')
        plt.xlabel('Read (bp)')

        # remove duplicates from legends
        handles, labels = currentAxis.get_legend_handles_labels()
        newLabels, newHandles = [], []
        for handle, label in zip(handles, labels):
            if label not in newLabels:
                newLabels.append(label)
                newHandles.append(handle)

        plt.legend(newHandles, newLabels)
        plt.savefig(output_plot)
    else:
        plt.savefig(output_plot)


if __name__ == '__main__':
    main()