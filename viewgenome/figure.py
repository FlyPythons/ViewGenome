from svgwrite import Drawing
from Bio import SeqIO

from .element import *


class Figure(object):

    def __init__(self):
        self.elements = []
        self.extra = []

        self.ticks = []

        self.svg = Drawing()

    def add(self, element, xy):

        self.elements.append(element)
        self.extra.append(xy)

        return element

    def add_genome(self, xy, name, length):

        genome = self.add(Genome(name, end=length), xy)

        return genome

    def add_genome_from_file(self, xy, file):

        file_format = file.split(".")[-1]

        if file_format in "gb|genbank":

            for record in SeqIO.parse(file, "genbank"):

                genome = self.add(Genome(name=record.id,
                                         end=len(record.seq),
                                         ), xy)

                for feature in record.features:
                    if feature.type == "CDS":
                        name = feature.qualifiers["locus_tag"][0]
                        genome.add_feature(name=name,
                                           start=feature.location.start,
                                           end=feature.location.end,
                                           strand=feature.location.strand,
                                           )

                break

            return genome

    def add_ticks(self, xy, ticks, labels):

        tick = Axis(ticks, labels)
        self.elements.append(tick)
        self.extra.append(xy)

        return tick

    def homology(self, name, color):

        homo = Homology(name, color)
        self.elements.append(homo)
        self.extra.append([0, 0])

        return homo

    def save(self, filename="out.svg"):

        scale = 1200.0 / max([i.length-i.break_length for i in self.elements if isinstance(i, Genome)])
        genomes = []

        for i, e in enumerate(self.elements):
            if isinstance(e, Genome):
                genomes.append(e.plot(scale=scale, x=self.extra[i][0], y=self.extra[i][1], stroke="black", stroke_width=2))
                continue

            self.svg.add(e.plot(scale=scale, x=self.extra[i][0], y=self.extra[i][1], stroke="black", stroke_width=2))

        for g in genomes:
            self.svg.add(g)

        self.svg.saveas(filename)


