#!/usr/bin/env python

from __future__ import print_function

import csv
import sys
import xml.etree.ElementTree as etree

import pyx
import pyx.connector

# Loading of the flux levels (one flux per line, one condition per column; first column is expected to be the gene name)
activity_fn = "PATH_TO/fluxes.csv"
activities, min_activity, max_activity = {}, sys.maxint, -sys.maxint
for entry in csv.reader(open(activity_fn, "rU")):
    entry_id, entry_activities = entry[0], [float(v) for v in entry[1:]]

    for timepoint_n, entry_activity in enumerate(entry_activities):
        activities.setdefault(timepoint_n, {})[entry_id] = entry_activity
        min_activity = min(min_activity, abs(entry_activity))
        max_activity = max(max_activity, abs(entry_activity))

# Loading of the Cytoscape XGMML export (version 1.1 is expected)
cytoscape_fn = "PATH_TO/pathway.xgmml"
cytoscape_g = etree.parse(cytoscape_fn)

cytoscape_ns = {
    "xgmml": "http://www.cs.rpi.edu/XGMML"}

graph_attrs = {}
graph_obj = cytoscape_g.getroot()

def get_attributes (node):
    return {entry.get("name"): entry.get("value") \
        for entry in node.findall("./xgmml:att", cytoscape_ns)}

# retrieve graph attributes
graph_attrs = get_attributes(graph_obj)

# list nodes
nodes = []
for node_obj in graph_obj.findall("./xgmml:node", cytoscape_ns):
    node_attrs = get_attributes(node_obj)
    graphics_obj = node_obj.find("./xgmml:graphics", cytoscape_ns)

    nodes.append((
        int(node_obj.get("id")),
        node_attrs["canonicalName"],
        float(graphics_obj.get("x")),
        -float(graphics_obj.get("y")),
        graphics_obj.get("type"),
        float(graphics_obj.get("w")),
        float(graphics_obj.get("h")),
        ))

# list edges
edges = []
for edge_obj in graph_obj.findall("./xgmml:edge", cytoscape_ns):
    edges.append((
        int(edge_obj.get("source")),
        int(edge_obj.get("target"))
        ))

#:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

pyx.unit.set(xscale = 50)

arrow_path = pyx.path.path(
    pyx.path.moveto(0, 0),
    pyx.path.lineto(2, 0),
    pyx.path.lineto(2, 0.5),
    pyx.path.lineto(3, -0.5),
    pyx.path.lineto(2, -1.5),
    pyx.path.lineto(2, -1),
    pyx.path.lineto(0, -1),
    pyx.path.closepath())

arrow_path = arrow_path.transformed(
    pyx.trafo.translate(-1.5, 0.5))

for timepoint_n in sorted(activities):
    canvas = pyx.canvas.canvas()

    # background shapes for nodes
    for (_, _, x, y, shape, w, h) in nodes:
        if (shape == "RECTANGLE"):
            shape_path = pyx.path.rect(x - w / 2, y - h / 2, w, h)
        elif (shape == "ELLIPSE"):
            shape_path = pyx.path.circle(x, y, w / 2)
        else:
            raise Exception("unknown shape: %s" % shape)

        canvas.draw(shape_path,
            [pyx.deco.filled(
                [pyx.color.cmyk.NavyBlue,
                 pyx.color.transparency(0.85)]),
             pyx.deco.stroked(
                [pyx.color.grey(0.8)])])

    for (node_id, node_label, x, y, _, _, _) in nodes:
        node_activity = activities[timepoint_n].get(node_label)
        if (node_activity is None):
            continue

        if (node_activity > 0):
            arrow_angle = 60
            arrow_color = pyx.color.rgb.green
        else:
            arrow_angle = -60
            arrow_color = pyx.color.rgb.red

        arrow_size = 10 + 50 * (abs(node_activity) - min_activity) / (max_activity - min_activity)

        canvas.draw(arrow_path,
            [pyx.trafo.scale(arrow_size),
             pyx.trafo.rotate(arrow_angle),
             pyx.trafo.translate(x, y),
             pyx.deco.filled(
                [arrow_color,
                 pyx.color.transparency(0.75)]),
             pyx.deco.stroked(
                [arrow_color,
                 pyx.color.transparency(0.15)])])

    # node labels
    node_to_path = {}
    for (node_id, node_label, x, y, _, _, _) in nodes:
        label_path = pyx.text.text(
            x, y, node_label,
            [pyx.text.halign.center,
             pyx.text.valign.middle])

        # canvas.fill(
        #     label_path.bbox().enlarged(2).path(),
        #     [pyx.color.rgb.white,
        #      pyx.color.transparency(0.25)])

        canvas.insert(label_path)
        node_to_path[node_id] = label_path

    # node connectors
    for (node_id_a, node_id_b) in edges:
        canvas.stroke(
            pyx.connector.line(
                node_to_path[node_id_a],
                node_to_path[node_id_b],
                boxdists = [10, 10]),
            [pyx.style.linewidth(2.5),
             pyx.deco.earrow(size = 12),
             pyx.color.grey(0.5)])

    canvas.writePDFfile("pathway_%d.pdf" % timepoint_n)