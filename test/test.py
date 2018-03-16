#!/usr/bin/env python
from viewgenome import figure


def main():

    fig = figure.Figure()

    # add genome to figure
    g1 = fig.add_genome_from_file([100, 150], "1.gb", )

    # add genome break to figure
    g1.add_break(5000, 10000)

    g2 = fig.add_genome_from_file([100, 300], "2.gb", )

    g3 = fig.add_genome_from_file([100, 450], "3.gb", )
    g3.add_break(10000, 15000)

    g4 = fig.add_genome_from_file([100, 600], "4.gb", )

    g5 = fig.add_genome_from_file([100, 750], "5.gb", )

    g6 = fig.add_genome(xy=[100, 900], name="Test", length=40000)

    # add feature to genome
    g6.add_feature(name="H1", start=1000, end=4000, strand=1)

    # add homology genes to figure

    h1 = fig.homology("h1", "blue")

    # use genome.feature[gene_name] to get gene
    h1.add_feature(g1.feature["wp2"],
                   g2.feature["PBC1_gp02"],
                   g3.feature["Fah02"],
                   g4.feature["CHERRY_0002"],
                   g5.feature["ISGA_2"])

    h1 = fig.homology(name="h2",
                      color="orange")
    h1.add_feature(g1.feature["wp41"],
                   g2.feature["PBC1_gp40"],
                   g3.feature["Fah35"],
                   g4.feature["CHERRY_0038"])

    # add ticks
    fig.add_ticks(xy=[100, 50],
                  ticks=[0, 2000, 10000],
                  labels=[0, 2000, 10000])
    fig.save()


if __name__ == "__main__":
    main()



