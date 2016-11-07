from imposm.parser import OSMParser
from copy import deepcopy
import math
import rtree
import argparse
import sys
import Dijkstra
import igraph 
import time
class OSMNode():
    def __init__(self, nodeID, lon, lat):
        self.nodeID = nodeID
        self.lon = lon
        self.lat = lat
    def __str__(self):
        return "OSMNode ID: %d \t Longitude: %f \t Latitude: %f\n" % (self.nodeID, self.lon, self.lat)


# A node is a OSMNode which lies in a intersection or is a dead end
class Node(OSMNode):
    def __init__(self, nodeID, lon, lat, interpolated=False):
        self.nodeID = nodeID
        self.lon = lon
        self.lat = lat
        self.interpolated = interpolated
        self.adjList = {} # adjacency list of the node with nodeID as key and the distance as value
        self.adjTimeList = {}  # adjacency list of the node with neighbor's nodeID as the key and the traversal time as value

    def getBB(self):
        return(self.lon, self.lat, self.lon, self.lat)

    def toString(self, flag):
        returnValue = ""
        for ids in self.adjList:
            if(flag == "d" or flag == "D"):
                returnValue += "%s \t %s \t %f\n" % (self.nodeID, ids, self.adjList[ids])
            elif(flag == "t" or flag == "T"):
                returnValue += "%s \t %s \t %s\n" % (self.nodeID, ids, str(self.adjTimeList[ids]))
            else:
                returnValue += "%s \t %s \t %f \t %s\n" % (self.nodeID, ids, self.adjList[ids], str(self.adjTimeList[ids]))
        return returnValue

    def __str__(self):
        return toString("")

# A class for making a graph
class Graph:
    def __init__(self, dbName=None):
        self.numVertices = 0
        self.numArcs = 0
        self.nodeDict = {}
        self.arcDict = {}
        self.duplicateArcDict = {}
        self.duplicateArcNum = 0
        if dbName==None:
            self.index = rtree.index.Index()
        else:
            self.index = rtree.index.Index(dbName)

    def addNode(self, node):
        if not self.nodeDict.has_key(node.nodeID):
            self.nodeDict[node.nodeID] = node
            self.numVertices += 1

    def addArc(self, srcID, dstID, dist, time):
        if not self.nodeDict.has_key(srcID):
            raise KeyError("The node with ID " + srcID + " is not present in the graph")
        if not self.nodeDict.has_key(dstID):
            raise KeyError("The node with ID " + dstID + " is not present in the graph")
        if not self.arcDict.has_key((srcID, dstID)):
            self.arcDict[(srcID, dstID)] = True
            self.nodeDict[srcID].adjList[dstID] = dist
            self.nodeDict[srcID].adjTimeList[dstID] = time
        else:
            if self.duplicateArcDict.has_key((srcID, dstID)):
                self.duplicateArcDict[(srcID, dstID)] += 1
            else:
                self.duplicateArcDict[(srcID, dstID)] = 1
            self.duplicateArcNum += 1

    def loadFromDB(self):
        for node in self.index.intersection((-180, -90, 180, 90), objects="raw"):
            print("entered first looop") 
            self.numVertices += 1
            self.nodeDict[node.nodeID] = node
            self.numArcs += len(node.adjList)
        for node in self.nodeDict:
            print("Entered second for loop")
            for edge in node.adjList:
                self.arcDict[(node.nodeID, edge)] = True
                self.addArc(node.nodeID, edge, node.adjList[edge], node.adjTimeList[edge])
                

    
    def buildRTree(self):
        for node in self.nodeDict.values():
            print("inside buildrtree")
            self.index.insert(long(node.nodeID), node.getBB(), obj=node)
    
    def toString(self, flag):
        returnValue = ""
        for node in self.nodeDict:
            returnValue += self.nodeDict[node].toString(flag)
        return returnValue

    def __str__(self):
        returnValue =""
        for node in self.nodeDict:
            thisNode = self.nodeDict[node]
            for ids in self.nodeDict[node].adjList:
                thatNode = self.nodeDict[ids]
                returnValue += "%s (%f, %f) \t %s (%f, %f) \t %f\n" %(thisNode.nodeID, thisNode.lon, thisNode.lat, thatNode.nodeID, thatNode.lon, thatNode.lat, thisNode.adjList[ids])

        return returnValue


class GraphDB:
    def __init__(self,dbName):
        self.index = rtree.index.Rtree(dbName)

    def getGraph(self, bb=(-180, -90, 180, 90)):
        ret = Graph()
        biggerBBox = (bb[0] - 0.01, bb[1] - 0.01 , bb[2] + 0.01 , bb[3] + 0.01)
        biggerGraph = Graph()

        for node in self.index.intersection(bb, 'raw'):
            ret.addNode(node)
            ret.numArcs += len(node.adjList)

        for node in self.index.intersection(biggerBBox, 'raw'):
            biggerGraph.addNode(node)
        interpolatedID =0
        i=0
        for node in ret.nodeDict.values():
            for edge in node.adjList.keys():
                if not ret.nodeDict.has_key(edge):
                    if biggerGraph.nodeDict.has_key(edge):
                        bbox = biggerGraph.nodeDict[edge].getBB()
                        lon = bbox[0]
                        lat = bbox[1]
                        x1 = node.getBB()[0]
                        y1 = node.getBB()[1]
                        slope = (lat - y1) / (lon - x1)
                        x2 =0
                        y2 =0
                        left = bb[0]
                        bottom = bb[1]
                        right = bb[2]
                        top = bb[3]
                        if (left <= lon) and (lon <= right):
                            if lat < bottom:
                                y2 = bottom
                            elif lat > top:
                                y2 = top
                            x2 = ((y2 -y1)/slope) + x1
                            if not ((left <= x2) and (x2 <= right) and (bottom <= y2) and (y2 <= top)): 
                                raise ValueError("Interpolation Error: 1")
                        else:
                            if (bottom <= lat) and (lat <= top):
                                if lon < left:
                                    x2 = left
                                elif lon > right:
                                    x2 = right
                                y2 = slope * (x2 -x1) + y1
                                if not ((left <= x2) and (x2 <= right) and (bottom <= y2) and (y2 <= top)):
                                    raise ValueError("Interpolation Error: 2")
                            else:
                                if (lon < left) and (lat > top):
                                    x2 = left
                                    y2 = slope * (x2 -x1) + y1
                                    if not ((bottom <= y2) and (y2 <= top)):
                                        y2 = top
                                        x2 = ((y2-y1)/slope) + x1
                                        if not((left <= x2) and (x2 <= top)):
                                            raise ValueError("Interpolation Error: 3")
                                elif (lon < left) and (lat < bottom):
                                    x2 = left
                                    y2 = slope * (x2 - x1) + y1
                                    if not ((bottom <= y2) and (y2 <= top)):
                                        y2 = bottom 
                                        x2 = ((y2-y1)/slope) + x1
                                        if not ((left <= x2) and (x2 <= top)):
                                            raise ValueError("Interpolation Error: 4")
                                elif ((lon > right) and (lat > top)):
                                    x2 = right
                                    y2 = slope * (x2- x1) + y1
                                    if not ((bottom <= y2) and (y2 <= top)):
                                        y2 = top
                                        x2 = ((y2-y1)/slope) + x1
                                        if not ((left <= x2) and (x2 <= top)):
                                            raise ValueError("Interpolation Error: 5")
                                elif ((lon > right) and (lat < bottom)):
                                    x2 = right
                                    y2 = slope * (x2- x1) + y1
                                    if not ((bottom <= y2) and (y2 <= top)):
                                        y2 = bottom
                                        x2 = ((y2-y1)/slope) + x1
                                        if not ((left <= x2) and (x2<= top)):
                                            raise ValueError("Interpolation Error: 6")
                    interpolatedNode = Node(str(i), x2, y2, True)
                    i+= 1
                    ret.addNode(interpolatedNode)
                    ret.addArc(node.nodeID, interpolatedNode.nodeID,-1, -1)
                    continue
                ret.addArc(node.nodeID, ret.nodeDict[edge].nodeID, node.adjList[edge], node.adjTimeList[edge])
        return ret
    ###







class Link():
    def __init__(self, linkID, refs, oneway, length, speedLimit):
        self.linkID = linkID
        #making deep copy of refs so that changes on refs does not affect the linkReferences
        self.linkReferences = deepcopy(refs)
        self.oneway = oneway
        self.length = length
        speedLimitParts = str(speedLimit).split(" ")
        # if the speedLimit contains only numbers then it is on kmph 
        if len(speedLimitParts)==1:
            # if there is not integer or float value for the speed limit then assign it to undefined
            if not (isinstance(speedLimit, float) or isinstance(speedLimit, int)):
                self.speedLimit = -1
            else:
                self.speedLimit = speedLimit
        # if the speedLimit contains both numeric value and unit then the conversion might be needed
        elif len(speedLimitParts)==2:
            if speedLimitParts[1] == 'km/h' or speedLimitParts[1]=='kmph':
                self.speedLimit = speedLimitParts[0]
            elif speedLimitParts[1] == 'mph':
                self.speedLimit = float(speedLimitParts[0]) * 1.60934
            elif speedLimitParts[1] == 'knots':
                self.speedLimit = float(speedLimitParts[0]) * 1.852
            else:
                self.speedLimit = -1
        else:
            self.speedLimit = -1

    # Calculating the distance between two nodes by using Haversine 
    #   formula and adding all the distances of adjacent nodes to 
    #   find out the total length of a link
    # The returned length is in kilometers
    @staticmethod
    def findLength(refs):
        i = 0
        length =0
        size = len(refs)
        # using while loop with index instead of for loop because we need to look at two nodes at a time
        while i < size -1:
            lon1 = math.radians(refs[i].lon)
            lon2 = math.radians(refs[i+1].lon)
            lat1 = math.radians(refs[i].lat)
            lat2 = math.radians(refs[i+1].lat)

            dlon = abs(lon2 - lon1)
            dlat = abs(lat2 - lat1)
            a = math.pow(math.sin(dlat/2),2) + math.cos(lat1) * math.cos(lat2)* math.pow(math.sin(dlon/2),2)
            c = 2 * math.atan2(math.sqrt(a), math.sqrt(1-a))
            #Using 6373 Km as radius of earth 
            d = 6373 * c
            length += d
            i += 1

        return length


    def __str__(self):
        retVal = '----------------\nLink ID: ' + str(self.linkID) +'\t Oneway: ' + self.oneway + '\t Length: ' + str(self.length) + '\nReferences: \n'
        strRef = ""
        for refs in self.linkReferences:
            strRef += str(refs)
        retVal += strRef + '\n------------\n'
        return retVal























class linkTracker():
    def __init__(self):
        self.highwayList = []
        self.referencesDict = {} 
        self.refinedLinkList = []
        self.referencesNodeDict = {}
        self.noSpeedLimitNum = 0


    # function for extracting highways from map
    # this function makes Link object for every highway and stores in highwayList
    #   and also stores all the node ids within the highways in recerecesDict as dictionary
    #   with node id as key and number of occurences of that node as value
    def extractHighway(self, ways):
        for osmid, tags, refs in ways:

            if 'highway' in tags:
                if (len(refs) > 1):
                    highwayID = (refs[0], refs[len(refs) -1])
                    refs = map(str, refs)
                    if tags['highway'] == 'motorway' or tags['highway'] == 'roundabout':
                        # if the highway has maxspeed tag then set speedLimit to the given speedLimit else set speedLimit to undefined
                        if 'maxspeed' in tags:
                            temp = Link(highwayID, refs, 'forward', 0, tags['maxspeed'])
                        else:
                            temp = Link(highwayID, refs, 'forward', 0, -1)
                    elif 'oneway' in tags:
                        # assigning the length of highways to zero as we
                        #   are not interested in their length
                        if tags['oneway'] == 'yes' or tags['oneway'] == 1 or tags['oneway'] == 'true':
                            if 'maxspeed' in tags:
                                temp = Link(highwayID, refs, 'forward', 0, tags['maxspeed'])
                            else: 
                                temp = Link(highwayID, refs, 'forward', 0, -1)
                        elif tags['oneway'] == 'no' or tags['oneway'] == 0 or tags['oneway'] == 'false' or tags['oneway'] == 'reversible':
                            if 'maxspeed' in tags:
                                temp = Link(highwayID, refs,'bidirectional',0, tags['maxspeed'])
                            else:
                                temp = Link(highwayID, refs, 'bidirectional', 0, -1)

                        else:
                            if 'maxspeed' in tags:
                                temp = Link(highwayID, refs, 'backwards', 0, tags['maxspeed'])
                            else:
                                temp = Link(highwayID, refs, 'backwards', 0, -1)
                    else:
                        if 'maxspeed' in tags:
                            temp = Link(highwayID, refs, 'bidirectional',0, tags['maxspeed'])
                        else:
                            temp = Link(highwayID, refs, 'bidirectional',0, -1)
                        
                    self.highwayList.append(temp)
                    for refID in refs:
                        if self.referencesDict.has_key(refID):
                            self.referencesDict[refID] += 1
                        else:
                            self.referencesDict[refID] = 1
        


    # this function gets coordinates for all the nodes associated with highways
    # This function will make node object for each node and store those node object in referencesNodeDict as dictionary
    #   with node id as key and the OSMNode object as value
    def getNodeData(self, coords):
        for osmid, lon, lat in coords:
            osmid = str(osmid)
            if osmid in self.referencesDict: 
                tempNode = OSMNode(str(osmid), lon, lat)
                self.referencesNodeDict[str(osmid)] = tempNode



    # this function breaks the highways in highwayList on the basis of intersections
    # and store them as link objects in refinedLinkList
    def refineLinks(self, graph, defaultSpeedLimit):
        # iterate through each highway in list
        for highway in self.highwayList:       
            # a temporary container for storing node objects 
            refs = []
            numRefs = len(highway.linkReferences)
            j =0
            # using while loop instead of for loop since we might not need to loop sequentially
            #   In the intersection we will have to consider the node in the intersection two times for two different links
            while(j < numRefs):
                id = highway.linkReferences[j]
                # append the OSMNode with given id in the temporary container
                refs.append(self.referencesNodeDict[id])

                # if the node id has appeared twice (meaning the node is intersection) or 
                #   the node id is the last reference of highway then create a new link with 
                #   whatever is in the temporary container of node id
                if self.referencesDict[id] > 1 or highway.linkReferences[len(highway.linkReferences) -1] == id:
                    # considering only the links which have more than one node as only one node represents intersection only
                    if len(refs) >1 :
                        finalIndex = len(refs) -1
                        linkID = (refs[0].nodeID, refs[finalIndex].nodeID)
                        linkLength = Link.findLength(refs)
                        # if the speedLimit of highway is -1 and defaultSpeedLimit is not -1 
                        #   then set speedLimit of highway to defaultSpeedLimit
                        if highway.speedLimit == -1 and defaultSpeedLimit != -1:
                            highway.speedLimit = defaultSpeedLimit
                        
                        #if the speedLimit is not undefined then calculate the traversal time else set traversal time to undefined
                        if not highway.speedLimit== -1 :
                            linkTraversalTime = linkLength/highway.speedLimit
                        else:
                            self.noSpeedLimitNum += 1
                            linkTraversalTime  = -1
                        temp = Link(linkID, refs, highway.oneway, linkLength, highway.speedLimit)
                        # first and last OSMNode of a link is a node because they are either intersection or dead end
                        # making node out of first OSMNode of the link and adding to the graph
                        node = Node(refs[0].nodeID, refs[0].lon, refs[0].lat)
                        graph.addNode(node)
                        

                        # making node out of last OSMNode of the link and adding to the graph
                        node2 = Node(refs[finalIndex].nodeID, refs[finalIndex].lon, refs[finalIndex].lat)
                        graph.addNode(node2)
                        # adding second node in the adjacent list of first node 
                        if(highway.oneway == 'forward' or highway.oneway == 'bidirectional'):
                            graph.addArc(refs[0].nodeID, refs[-1].nodeID, linkLength, linkTraversalTime)
                        if(highway.oneway == 'backwards' or highway.oneway == 'bidirectional'):
                            # adding first node in the adjacent list of the second node
                            graph.addArc(refs[-1].nodeID, refs[0].nodeID, linkLength, linkTraversalTime)
                        self.refinedLinkList.append(temp)
                        del refs[:]
                        
                        # if the node is in intersection then continue without increasing the value of j
                        # i.e start with the same node for next link
                        if self.referencesDict[id] > 1:
                            continue
                j+= 1        
    
def plot(graph, program_start_time):
    startTime = time.time()
    iGraph = igraph.Graph(directed=False)

    iGraph.add_vertices(len(graph.nodeDict))
    iGraph.vs["name"] = graph.nodeDict.keys()
    iGraph.add_edges(graph.arcDict.keys())

    mylayout = []

    for keys in graph.nodeDict:
        mylayout.append((graph.nodeDict[keys].lon, -(graph.nodeDict[keys].lat)))

   

    visual_style = {}
    visual_style["vertex_size"] =3
    visual_style["vertex_color"] = "blue"
    visual_style["edge_color"] = "red"
    visual_style["edge_curved"] = 0
    visual_style["layout"] = mylayout
    endTime = time.time()
    dt = endTime - startTime
    dtProgramTime = endTime - program_start_time
    sys.stderr.write("Percentage of time required to build igraph: " + str((dt/dtProgramTime)*100)+ "%\n")

    igraph.plot(iGraph, **visual_style)
    return endTime



def main():
    programStartTime = time.time()
    argparser = argparse.ArgumentParser()
    argparser.add_argument("pbf_file", help="The pbf file that needs to be converted")
    argparser.add_argument("--distance", "-d", help="Outputs the distance between two nodes", action="store_true")
    argparser.add_argument("--traversal_time" , "-t", help="Outputs the traversal time between two nodes", action="store_true")
    argparser.add_argument("--node", "-n", help="Outputs the information about nodes only. <NodeID>    <Latitude>   <Longitude>", action="store_true")
    argparser.add_argument("-l", dest="defaultSpeedLimit",  type=float ,  help="Sets default speed limit(in Kmph) for segments without speed limit tag", default=-1)
    argparser.add_argument("--output", "-o", type=str, help="Specifies the output file. Default to terminal", default="")
    
    args = argparser.parse_args()
    
    if ((args.defaultSpeedLimit) != -1 and (args.defaultSpeedLimit <= 0)):
        raise argparse.ArgumentTypeError("Default Speed Limit  not in range.(DefaultSpeedLimit > 0)")
    
    tracker = linkTracker()
    extractHighwayStart = time.time()
    p = OSMParser(ways_callback=tracker.extractHighway)
    print("Collecting information about highways")
    #calling parse to collect information about highways
    p.parse(args.pbf_file)
    extractHighwayEnd = time.time()
    dtExtractHighway = extractHighwayEnd - extractHighwayStart
    
    getNodeDataStart = time.time()
    q = OSMParser(coords_callback=tracker.getNodeData)
    print("Collecting nodes informaion")
    #calling parse to collect information about the nodes associated with highways
    q.parse(args.pbf_file)
    getNodeDataEnd = time.time()
    dtNodeData = getNodeDataEnd - getNodeDataStart

    graph = Graph("uno") 

    #plot(graph, layout=layout)

    print("Refining highways")
    refineHighwayStart = time.time()
    #calling refineLinks to break down highways on the basis of intersections
    tracker.refineLinks(graph, args.defaultSpeedLimit)
    refineHighwayEnd = time.time()
    dtRefineHighway = refineHighwayEnd - refineHighwayStart

    programEndTime = plot(graph, programStartTime)
    dtProgramTime = programEndTime - programStartTime
 
    if(args.output != ""):
        sys.stdout = open(args.output, 'w')
   

    graph.buildRTree()
    if(args.node):
        for id in graph.nodedict:
            print("%d \t %f \t %f"%(id, graph.nodedict[id].lat, graph.nodedict[id].lon))
    elif(args.distance and (not args.traversal_time)):
        print(graph.tostring("d"))
    elif((not args.distance) and args.traversal_time):
        print(graph.tostring("t"))
    else:
        print(graph.toString(""))

    
    #dj = Dijkstra.Dijkstra(graph, 115343360)
    #print("Applying Dijkstra's Alogrithm")
    #dj.dijkstra()
    #value = dj.getValue(10)
    #print(dj.getValue(10))
    sys.stderr.write("Percentage of time taken to extract highway information: "+str((dtExtractHighway/dtProgramTime)*100) + "%\n")
    sys.stderr.write("Percentage of time taken to extract node data: " + str((dtNodeData/dtProgramTime)*100) + "%\n")
#   sys.stderr.write("Percentage of time taken to refine highway: " + str((dtRefineHighway/dtProgramTime)*100) + "%\n")

    sys.stderr.write("Total Links: " + str(len(graph.arcDict))+ "\n")
    sys.stderr.write("Total OSM nodes: " + str(len(tracker.referencesNodeDict)) + "\n")
    sys.stderr.write("Total Nodes/Vertices: " + str(graph.numVertices) + "\n")
    sys.stderr.write("Total links without speedLimit tag: " + str(tracker.noSpeedLimitNum) + "\n")
    sys.stderr.write("Total number of duplicate arcs: " + str(graph.duplicateArcNum) + "\n")
    sys.stderr.write("Number of distinct duplicate arcs: " + str(len(graph.duplicateArcDict)) + "\n")
   #print("Total Links: " + str(len(tracker.refinedLinkList)))
   #print("Total node objects: " + str(len(tracker.referencesNodeDict)))
   #print("Total vertices in graph: " + str(graph.numVertices))   

if __name__ == '__main__': 
    main()
