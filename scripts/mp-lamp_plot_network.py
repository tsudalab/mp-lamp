#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import math
import networkx as nx
import pandas as pd
import sys, re

CHR2CORLR_DIC = {
    "1": "#66C2A5", # 色はbrewer.pal(5, "Set2")を利用
    "2": "#FC8D62",
    "3": "#8DA0CB",
    "4": "#E78AC3",
    "5": "#A6D854"
}

def buildNetwork(file):
    #data = pd.read_table(file, sep="\s+")
    data = readMPLAMPoutput(file)

    G = nx.Graph()
    for vitems in data.COMB:
        for i in range(len(vitems)):
            if G.has_node(vitems[i]):
                G.node[vitems[i]]["count"] += 1
            else:
                G.add_node(vitems[i], count = 1)
                #G.add_node(vitems[i], {"count":1})
        for i in range(len(vitems)):
            for j in range(i+1, len(vitems)):
                if G.has_edge(vitems[i], vitems[j]):
                    G[vitems[i]][vitems[j]]["weight"] += 1
                    #G.edge[vitems[i]][vitems[j]]["weight"] += 1
                else:
                    G.add_edge(vitems[i], vitems[j], weight = 1)
                    #G.add_edge(vitems[i], vitems[j], {"weight":1})

    maxcount = 1
    if len(G.nodes()) > 0:
        maxcount = max([G.node[node]['count'] for node in G.nodes()])
    #maxcount = max([G.node[node]['count'] for node in G.nodes()])
    maxweight = 1
    if len(G.edges()) > 0:
        maxweight = max([G[edge[0]][edge[1]]['weight'] for edge in G.edges()])
    #maxweight = max([G.edge[edge[0]][edge[1]]['weight'] for edge in G.edges()])
    for node in G.nodes():
        G.node[node]['count'] /= float(maxcount)
        #G.node[node]['count'] /= maxcount
    for edge in G.edges():
        G[edge[0]][edge[1]]['weight'] /= float(maxweight)
        #G.edge[edge[0]][edge[1]]['weight'] /= maxweight

    return(G)

def plotNetwork(nx_graph,
                output = None,
                fig_width = 1000,
                dpi = 72,
                k = None,
                node_size = None,
                edge_width = None,
                fixed_node_size = True,
                rs2chr_dict = {}):
    fig_width = fig_width / 72
    plt.figure(figsize = (fig_width, fig_width), dpi = 72)
    pos = nx.spring_layout(nx_graph, k = k)

    if node_size is None:
        node_size = fig_width * fig_width * 72 * 72 / 20 / len(nx_graph.nodes())
    if edge_width is None:
        edge_width = 2 * math.sqrt(node_size / math.pi)
        if len(nx_graph.edges()) > 0:
            edge_width = 2 * math.sqrt(node_size / math.pi) / math.sqrt(len(nx_graph.edges()))
        #edge_width = 2 * math.sqrt(node_size / math.pi) / math.sqrt(len(nx_graph.edges()))
    if not fixed_node_size:
        node_size = [d['count'] * node_size for (n,d) in nx_graph.nodes(data = True)]
    edge_width = [d['weight'] * edge_width for (u,v,d) in nx_graph.edges(data = True)]

    # Node color
    ncolor_list = []
    for n, d in nx_graph.nodes(data = True):
        if n in rs2chr_dict.keys():
            chrid = rs2chr_dict[n]
            ncolor_list.append(CHR2CORLR_DIC[chrid])
        else:
            ncolor_list.append('lightgrey')

    nx.draw_networkx_edges(G, pos, alpha=0.4, edge_color='grey', width=edge_width)
    nx.draw_networkx_nodes(G, pos, node_color=ncolor_list, alpha=.9, node_size=node_size)
    nx.draw_networkx_labels(G, pos, fontsize=14, font_weight="bold", alpha=0.5)

    plt.axis('off')
    if output is None:
        plt.show()
    else:
        plt.savefig(output, dpi = dpi)
    plt.close()


def readMapFile(map_filename):
    """
    Read a MAP file,
    return the dictionary whose key and values are SNP ID and chromosome number, respectively.
    """
    fi = open(map_filename, 'r')
    rs2chr_dict = {}
    for line in fi:
        s = re.split('\s+', line.strip())
        rsid = s[1]
        chrid = s[0]
        rs2chr_dict[rsid] = chrid
    fi.close()
    return rs2chr_dict

def readMPLAMPoutput(mplamp_filename):
    """
    Read mp-lamp output file
    """
    temp = []

    for line in open(mplamp_filename):
        li = line.strip()
        if not li.startswith("#"):
            temp.append(li)

    df = pd.DataFrame(temp, columns=['RAW_PVAL'])

    cols = ['RAW_PVAL','CORR_PVAL','FREQ','POS','NU_ITEM','COMB']
    df[cols] = df.RAW_PVAL.str.split('\t', n=5, expand=True)
    df.COMB = df.COMB.str.split('\t')
    return df

if __name__ == "__main__":
    import argparse
    import os

    parser = argparse.ArgumentParser(
        description='Draw Network graph from LAMP result...')
    parser.add_argument('-v', '--version',
                        action='version',
                        version='%(prog)s 0.9.0')
    parser.add_argument('lamp_result',
                        action='store',
                        nargs=None,
                        const=None,
                        default=None,
                        type=str,
                        choices=None,
                        help='LAMP result file.',
                        metavar=None)
    parser.add_argument('output_image',
                        action='store',
                        nargs=None,
                        const=None,
                        default=None,
                        type=str,
                        choices=None,
                        help='Output image file.',
                        metavar=None)
    parser.add_argument('--fig-size',
                        action='store',
                        nargs=None,
                        const=None,
                        default=1000,
                        type=int,
                        choices=None,
                        help='image size.(pixels, default:1000)',
                        metavar='INT')
    parser.add_argument('--dpi',
                        action='store',
                        nargs=None,
                        const=None,
                        default=72,
                        type=int,
                        choices=None,
                        help='image dpi.(default:72)',
                        metavar='INT')
    parser.add_argument('-k',
                        action='store',
                        nargs=None,
                        const=None,
                        default=None,
                        type=float,
                        choices=None,
                        help='optimal distance between nodes.(default:None)',
                        metavar='FLOAT')
    parser.add_argument('--node-size',
                        action='store',
                        nargs=None,
                        const=None,
                        default=None,
                        type=int,
                        choices=None,
                        help='maximum node area size.',
                        metavar='INT')
    parser.add_argument('--edge-width',
                        action='store',
                        nargs=None,
                        const=None,
                        default=None,
                        type=int,
                        choices=None,
                        help='maximum edge width.',
                        metavar='INT')
    parser.add_argument('--node-size-fixed',
                        action='store_true',
                        help='plot fixed node size.')
    parser.add_argument('--rep',
                        action='store',
                        nargs=None,
                        const=None,
                        default=1,
                        type=int,
                        choices=None,
                        help='save image repeatedly.',
                        metavar='INT')
    parser.add_argument('--map-file',
                        action='store',
                        nargs=None,
                        const=None,
                        default=None,
                        type=str,
                        choices=None,
                        help='MAP filename')

    args = parser.parse_args()

    lampResultFile = args.lamp_result
    outputFile = args.output_image
    figSize = args.fig_size
    dpi = args.dpi
    nodeSize = args.node_size
    edgeWidth = args.edge_width
    k = args.k
    fixed = args.node_size_fixed
    n = args.rep
    map_filename = args.map_file

    # read the MAP file to color node according to Chr
    rs2chr_dict = {}
    if not (map_filename == None):
        rs2chr_dict = readMapFile(map_filename)

    G = buildNetwork(lampResultFile)
    # If the graph is empty,
    # output the error message and stop this script
    if len(G.nodes()) == 0:
        sys.stderr.write("%s does not contain a node\n" % lampResultFile)
        sys.exit()

    (name, ext) = os.path.splitext(outputFile)
    if ext.lower() not in ('.png', '.pdf', '.ps', '.eps', '.svg'):
        name = outputFile
        ext = ".png"
    if n > 1:
        path_list = [''.join(map(str, (name, '_', i+1, ext))) for i in range(n)]
    else:
        path_list = [''.join(map(str, (name, ext)))]

    for out in path_list:
        plotNetwork(G, output = out, fig_width = figSize, dpi = dpi, node_size = nodeSize, k = k, fixed_node_size = fixed, rs2chr_dict = rs2chr_dict)
