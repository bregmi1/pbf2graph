from Map2Graph import *

def plot(graph):
    startTime = time.time()
    iGraph = igraph.Graph(directed=False)

    iGraph.add_vertices(len(graph.nodeDict))
    iGraph.vs["name"] = graph.nodeDict.keys()
    iGraph.add_edges(graph.arcDict.keys())

    mylayout = []
    interpolated = []

    for keys in graph.nodeDict:
        mylayout.append((graph.nodeDict[keys].lon, -(graph.nodeDict[keys].lat)))
        interpolated.append((graph.nodeDict[keys].interpolated))
    
    iGraph.vs["interpolated"] = interpolated
   

    visual_style = {}
    visual_style["vertex_size"] =4
    color_dict = {True:"red", False:"blue"}

    visual_style["vertex_color"] = [color_dict[value] for value in iGraph.vs["interpolated"]]
    visual_style["edge_color"] = "red"
    visual_style["edge_curved"] = 0
#   visual_style["vertex_label"] = [long(id)%100 for id in iGraph.vs["name"]]
    visual_style["vertex_label_size"] = 10
    visual_style["vertex_label_dist"] = 2
    visual_style["layout"] = mylayout
    visual_style["target"] = "output.png"
    endTime = time.time()
    dt = endTime - startTime
    sys.stderr.write("time required to build igraph: " + str(dt)+ " seconds\n")
    igraph.plot(iGraph, **visual_style)

def main():
    bbox = (-90.070343,30.023075, -90.060620, 30.032894)
    graphDB = GraphDB("nola")
    graph = graphDB.getGraph(bbox)
    plot( graph)


if __name__=='__main__':
    main()
