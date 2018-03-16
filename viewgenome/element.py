from svgwrite import Drawing


class Genome(object):

    element = "Genome"

    def __init__(self, name, start=0, end=0, **extra):
        assert isinstance(start, int)
        assert isinstance(end, int)
        assert end > start

        self.name = name
        self.start = start
        self.end = end
        self.length = end - start

        self._features = {}

        self.position = None
        self._break = []

        self.extra = {"stroke": "black",
                      "stroke_width": 2
                      }
        self.extra.update(extra)

    @property
    def features(self):

        return sorted(self._features.values(), key=lambda v: v.start)

    @property
    def feature(self):
        return self._features

    @property
    def break_length(self):
        r = 0
        for i in self._break:
            r += i[1] - i[0]
        return r

    def add_feature(self, name, start, end, strand):

        assert name not in self._features

        feature = Feature(name, start, end, strand)

        self._features[name] = feature
        feature.parent[self.element] = self

        return feature

    def add_break(self, start, end):

        assert isinstance(start, int)
        assert isinstance(end, int)
        assert start < end

        self._break.append((start, end))

    def get_fragments(self):

        r = []
        prev_break = [0, 0]
        i = 0

        for i in range(len(self._break)):

            b = self._break[i]
            frag_genome = Genome(name="%s-%s" % (self.name, i),
                                 start=prev_break[1],
                                 end=b[0])

            prev_break = b
            r.append(frag_genome)

        frag_genome = Genome(name="%s-%s" % (self.name, i+1),
                             start=prev_break[1],
                             end=self.length)

        r.append(frag_genome)

        for k in self._features.keys():
            g = self._features[k]

            for f in r:

                if g.start >= f.start and g.end <= f.end:
                    self._features[k].frag_genome = f
                    break

        return r

    def plot(self, x, y, scale, break_width=10, **extra):

        svg = Drawing()

        break_length = 0
        fragments = self.get_fragments()
        g = svg.g()

        for i in range(len(fragments)):

            f = fragments[i]

            if i > 0:
                break_length += f.start - fragments[i-1].end - break_width/scale

            p1 = [x + (f.start - break_length)*scale, y]
            p2 = [x + (f.end - break_length)*scale, y]

            f.position = p1

            line1 = svg.line(start=p1, end=p2, **extra)

            line2 = svg.line(start=[p1[0], p1[1] - 5], end=[p1[0], p1[1] + 5], **extra)
            text2 = svg.text(text=f.start,
                             insert=(p1[0], p1[1] + 10),
                             font_size=8, font_family="Arial",
                             transform='rotate(%s, %s, %s)' % (90, p1[0], p1[1] + 10))
            line3 = svg.line(start=[p2[0], p1[1] - 5], end=[p2[0], p1[1] + 5], **extra)
            text3 = svg.text(text=f.end,
                             insert=(p2[0], p2[1] + 10),
                             font_size=8, font_family="Arial",
                             transform='rotate(%s, %s, %s)' % (90, p2[0], p1[1] + 10))
            g.add(line1)
            g.add(line2)
            g.add(text2)
            g.add(line3)
            g.add(text3)

        text1 = svg.text(
            text=self.name, insert=(10, y),
            font_size=10, font_family="Arial"
        )
        g.add(text1)
        text2 = svg.text(
            text="%s nt" % self.length,
            insert=(x / 2.0 - 10, y + 12),
            font_size=8, font_family="Arial"
        )
        g.add(text2)

        for e in self.features:
            if e.frag_genome:
                patch, text = e.plot(scale=scale)
                g.add(patch)
                g.add(text)

        return g

    def __repr__(self):
        return "{0.name!s}-{0.length!s}".format(self)


class Feature(object):

    def __init__(self, name, start, end, strand, **extra):

        self.name = name
        self.start = start
        self.end = end
        self.strand = strand
        self.length = end - start + 1

        self.parent = {"Genome": None,
                       "Homology": None
                       }

        self.frag_genome = None
        self.extra = extra
        self.position = None
        self.height = 15

    def plot(self, shape="arrow", strand=True, scale=1/2000, **extra):

        height = self.height

        if strand:
            y1 = self.frag_genome.position[1] - height / 2.0 - self.strand * height
            y2 = y1 + height
        else:
            y1 = self.frag_genome.position[1] - height / 2.0
            y2 = self.frag_genome.position[1] + height / 2.0

        if shape == "rect":
            p1 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale, y1]
            p2 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale, y1]
            p3 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale, (y1 + y2) / 2]
            p4 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale, y2]
            p5 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale, y2]

        if shape == "arrow":

            if self.strand == 1:

                p1 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale, y1]
                p2 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale - height/3, y1]
                p3 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale, (y1 + y2)/2]
                p4 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale - height/3, y2]
                p5 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale, y2]
            else:
                p1 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale, y1]
                p2 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale + height/3, y1]
                p3 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale, (y1 + y2)/2]
                p4 = [self.frag_genome.position[0] + (self.start - self.frag_genome.start)*scale + height/3, y2]
                p5 = [self.frag_genome.position[0] + (self.end - self.frag_genome.start)*scale, y2]

            if self.strand * p2[0] < self.strand * p1[0]:
                p2 = p1
                p4 = p5

        self.position = p1
        if self.strand == -1:
            self.position = p3

        svg = Drawing()

        if self.parent["Homology"]:
            patch = svg.polygon(
                points=[p1, p2, p3, p4, p5],
                fill=self.parent["Homology"].color
            )
        else:
            patch = svg.polygon(
                points=[p1, p2, p3, p4, p5],
                fill="white", stroke="black"
            )

        t = [(p1[0] + p3[0]) / 2, p1[1] - 1]
        a = 1

        if strand and self.strand == -1:
            t = [(p1[0] + p3[0]) / 2, p5[1] + 5]
            a = -1

        text = svg.text(text=self.name, insert=t, font_size="8",
                        transform='rotate(%s, %s, %s)' % (-45 * a, t[0], t[1]))

        return patch, text

    def __repr__(self):
        return "{0.parent!s}-{0.name!s}".format(self)


class Homology(object):

    element = "Homology"

    def __init__(self, name, color="black"):

        self.name = name
        self.features = []
        self.color = color

    def add_feature(self, *feature):

        for f in feature:
            self.features.append(f)
            f.parent[self.element] = self

    def sort(self, *order):

        return sorted(self.features, key=lambda d: order.index(d.genome))

    def plot(self, shape="rect", scale=0, strand=True, **extra):

        svg = Drawing()
        g = svg.g()

        previous = None

        for feature in self.features:

            if not previous:
                previous = feature
                continue

            x1, y1 = previous.position
            x2, y2 = feature.position

            y1 += (2 + previous.strand)*previous.height + previous.height/2

            y2 += (feature.strand - 1) * feature.height - previous.height / 2

            p1 = [x1, y1]
            p2 = [x1 + previous.length*scale, y1]
            p3 = [x2 + feature.length*scale, y2]
            p4 = [x2, y2]

            if shape == "rect":
                g.add(svg.polygon(
                    points=[p1, p2, p3, p4],
                    fill=self.color,
                    fill_opacity=0.6
                ))

            if shape == "line":
                g.add(svg.polyline(
                    [((p1[0]+p2[0])/2, y1-previous.height/2),
                     ((p1[0] + p2[0]) / 2, y1),
                     ((p3[0]+p4[0])/2, y2),
                     ((p3[0] + p4[0]) / 2, y2+feature.height/2),
                     ],
                    fill_opacity=0,
                    stroke=self.color,
                    stroke_width=2
                ))

            previous = feature

        return g


class Axis(object):

    element = "Axis"

    def __init__(self, ticks, labels):
        assert len(ticks) == len(labels)

        self.ticks = ticks
        self.labels = labels

    def plot(self, x, y, scale, **extra):

        assert isinstance(scale, (int, float))
        assert isinstance(x, (int, float))
        assert isinstance(y, (int, float))

        svg = Drawing()
        g = svg.g()

        p1 = [x, y]
        p2 = [x + max(self.ticks) * scale, y]

        g.add(svg.line(start=p1, end=p2,
                       stroke="black", stroke_width=1
                       ))

        for i, tick in enumerate(self.ticks):
            p1 = [x + tick * scale, y - 5]
            p2 = [x + tick * scale, y + 0.5]

            g.add(svg.line(start=p1, end=p2,
                           stroke="black", stroke_width=1
                           ))

            g.add(svg.text(text="%s" % self.labels[i],
                           insert=[x + tick * scale, y+10],
                           font_size=8, font_family="Arial"))

        return g


