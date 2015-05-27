import numpy
from matplotlib import pyplot

def plot_embedding(points, labels, colors=None):
    points = numpy.array(points)
    figure = pyplot.figure()
    pyplot.axes([0., 0., 1., 1.])
    pyplot.scatter(points[:, 0], points[:, 1], color=colors)
    for i in range(len(labels)):
        pyplot.annotate(
            labels[i], points[i], fontsize=4, horizontalalignment='left')
    return [figure]

