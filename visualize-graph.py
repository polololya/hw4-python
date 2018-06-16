from Bio import SeqIO
from graphviz import Digraph
import argparse
import pydot

class Vertex:

    def __init__(self, seq):
        self.seq = seq
        self.vertex_coverage = 1
        self.in_edges = {}
        self.out_edges = {}

    def increase_vertex_coverage(self):
        self.vertex_coverage += 1

class Edge:

    def __init__(self, current_kmer, next_kmer):
        self.seq = current_kmer + next_kmer[-1]
        self.in_vertices = {}
        self.out_vertices = {}

        self.edge_coverage = 0

    def increase_edge_coverage(self):
        self.edge_coverage += 1

    def calculation_edge_coverage(self, prev_vertex_cov, next_vertext_cov):
        self.edge_coverage = (prev_vertex_cov + next_vertext_cov)/2

class Graph:

    def __init__(self, k):
        self.vertices = {}
        self.k = k

    def add_read(self, read):

        if len(read) < self.k:
            return

        # first k-mer
        kmer = read[:self.k]
        if kmer in self.vertices:
            self.vertices[kmer].increase_vertex_coverage()
        else:
            self.vertices[kmer] = Vertex(kmer)

        # next k-mer iterations:
        for i in range(1, len(read)-self.k+1, 1):
            next_kmer = read[i:i+self.k]
            if next_kmer in self.vertices:
                self.vertices[next_kmer].increase_vertex_coverage()
            else:
                self.vertices[next_kmer] = Vertex(next_kmer)

            # add new edge
            new_edge = Edge(kmer, next_kmer)
            # add vertices
            self.vertices[next_kmer].in_edges[kmer] = [new_edge]
            self.vertices[kmer].out_edges[next_kmer] = [new_edge]
            kmer = next_kmer

    def coverage_calculating(self):
        for current_vertex in self.vertices.keys():
            for next_vertex in self.vertices[current_vertex].out_edges.keys():
                self.vertices[current_vertex].out_edges[next_vertex][0]\
                      .calculation_edge_coverage(self.vertices[current_vertex].vertex_coverage,
                      self.vertices[next_vertex].vertex_coverage)

    def visualize_graph(self, arg):
        dot = Digraph(comment='Assemble')
        if arg == 'full':
            for k, v in self.vertices.items():
                dot.node(k, label='{}'.format(k))
                for kk, vv in v.out_edges.items():
                    dot.edge(k, kk, label='{}'.format(vv[0].seq))
        else:
            for k, v in self.vertices.items():
                dot.node(k, label='cov={}'.format(v.vertex_coverage))
                for kk, vv in v.out_edges.items():
                    dot.edge(k, kk, label='coverage:{} length:{}'.format(vv[0].edge_coverage,len(vv[0].seq)))
        print(dot.source)
        dot.view()
        dot.save()





if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', '--input',type=str)
    parser.add_argument('-k', '--kmer', help='k-mer length', default=15, type=int)
    parser.add_argument('-t', '--type', help='visualization type: full/short', default='full', type=str)
    parser.add_argument('-s', '--strand', help='fw|bw', default='fw', type=str)

    args = parser.parse_args()
    my_graph = Graph(args.kmer)
    if args.strand == 'fw':
        with open(args.input) as f:
            for record in SeqIO.parse(f, 'fastq'):
                my_graph.add_read(str(record.seq))
    else:
        with open(args.input) as f:
            for record in SeqIO.parse(f, 'fastq'):
                my_graph.add_read(str(record.reverse_complement().seq))

    my_graph.coverage_calculating()
    my_graph.visualize_graph(arg=args.type)