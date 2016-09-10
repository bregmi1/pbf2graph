from random import randint

class Dijkstra:
    def __init__(self, graph, initial):
        self.initial = initial
        self.visited = {self.initial:0}
        self.graph = graph
        self.path = {}
        

    def dijkstra(self):
        nodes = set()
        for nodeID in self.graph.nodeDict:
            nodes.add(nodeID)

        while nodes:
            min_node = None
            
            for node in nodes:
                if node in self.visited:
                    if min_node is None:
                        min_node = node
                    elif self.visited[node] < self.visited[min_node]:
                        min_node = node
            if min_node is None:
                break

            nodes.remove(min_node)
            current_weight = self.visited[min_node]


            for edge in self.graph.nodeDict[min_node].adjList:
                weight = current_weight + self.graph.nodeDict[min_node].adjList[edge]
                if edge not in self.visited or weight < self.visited[edge]:
                    self.visited[edge] = weight
                    self.path[edge] = min_node
    
    def getValue(self, id):
        self.__str__(); 
    
    def __str__(self):
        for ids in self.visited:
            print("%s \t %s \t %f" %(self.initial, ids, self.visited[ids]))



class Node:
    def __init__(self, nodeID):
        self.nodeID = nodeID
        self.adjList = {}

class Graph:
    def __init__(self):
        self.nodeDict = {}
        self.arcDict = {}

    def addNode(self, node):
        if not self.nodeDict.has_key(node.nodeID):
            self.nodeDict[node.nodeID] = node

    def addArc(self, srcNode, dstNode, distance):
        if not self.arcDict.has_key((srcNode.nodeID, dstNode.nodeID)):
            self.arcDict[(srcNode.nodeID, dstNode.nodeID)] = True
            self.nodeDict[srcNode.nodeID].adjList[dstNode.nodeID] = distance


def main():

    a = Node("a")
    b = Node("b")
    c = Node("c")
    d = Node("d")
    e = Node("e")


    gr = Graph()
    gr.addNode(a)
    gr.addNode(b)
    gr.addNode(c)
    gr.addNode(d)
    gr.addNode(e)


    gr.addArc(a,b,5)
    gr.addArc(a, c, 20)
    gr.addArc(b, d, 6)
    gr.addArc(b, e, 4)
    gr.addArc(d, e, 2)
    gr.addArc(e, c, 3)
    gr.addArc(e, d, 1)

    dj = Dijkstra(gr, "a")
    dj.dijkstra()
    #print(dj.__str__())
    dj.getValue(23)

if __name__ == '__main__':
    main()
