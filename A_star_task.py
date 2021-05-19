import re
import copy

flatten = lambda t: [item for sublist in t for item in sublist]

def read_file(path):
    """
    Reads file at path, the output will be a dictionary with:
    {   data: the matrix of data,
        upset: the list of pairs of upset kids
        questioned: the list of quoestioned kids
        message: the starting and goal point of message}
    """
    f = open(path, "r")
    d = f.read()
    data = {}
    data["data"], d = re.split("suparati\n|ascultati:\n|mesaj: ", d, 1)
    data["upset"], d = re.split("suparati\n|ascultati:\n|mesaj: ", d, 1)
    data["questioned"], data["message"] = re.split("suparati\n|ascultati:\n|mesaj: ", d, 1)
    for key in data:
        data[key] = data[key].split("\n")
        if not data[key][-1]:
            data[key].pop()
    data["message"][0] = data["message"][0].split(" -> ")
    data["message"] = flatten(data["message"])
    for key in data:
        for j in range(len(data[key])):
            data[key][j] = data[key][j].split(" ")
    data["questioned"] = flatten(data["questioned"])
    data["message"] = flatten(data["message"])
    return data

read_file("./input/1_in.in")

class NodeTransition:
    def __init__(self, index, id, currentQuestioned, position, father, cost):
        self.index = index
        self.id = id
        self.currentQuestioned = currentQuestioned
        self.position = position
        self.father = father
        self.g = cost

    def getPath(self):
        path = [self.id]
        node = self.father
        while node is not None:
            path.insert(0, node.id)
            node = node.father
        return path

    def showPath(self):
        strPath = str(self.id)
        node = self.father
        prev = self
        while node is not None:
            if prev.position["y"] < node.position["y"]:
                if prev.position["y"] // 2 == node.position["y"] // 2:
                    strPath = node.id + " < " + strPath
                elif prev.position["y"] // 2 < node.position["y"] // 2:
                    strPath = node.id + " << " + strPath
            elif prev.position["y"] > node.position["y"]:
                if prev.position["y"] // 2 == node.position["y"] // 2:
                    strPath = node.id + " > " + strPath
                elif prev.position["y"] // 2 > node.position["y"] // 2:
                    strPath = node.id + " >> " + strPath
            elif prev.position["x"] < node.position["x"]:
                strPath = node.id + " ^ " + strPath
            else: 
                strPath = node.id + " v " + strPath
            prev = node
            node = node.father
        return strPath

    def isInPath(self, nodeId):
        path = self.getPath()
        if nodeId in path:
            return True
        return False

    def __repr__(self):
        myStr = ""
        myStr += self.id + "(path = "
        myStr += self.showPath()
        myStr += ", cost = " + str(self.g)
        return myStr



class Graph:
    movesX = [0, 0, 1, -1]
    movesY = [1, -1, 0, 0]
    def __init__(self, nodes, matrix, start, goal, questioned, upset, timeQuestioned):
        self.nodes = nodes
        self.matrix = matrix
        self.start = start
        self.goal = goal
        self.questioned = questioned
        self.upset = upset
        self.timeQuestioned = timeQuestioned

    def getIndexNode(self, name):
        """
        self.nodes will have the format node: location {x: , y: }
        since every node-name is unique
        """
        return self.nodes[name]

    def test_goal(self, currentNode):
        return self.goal == currentNode.id

    def generateSuccessors(self, currentNode):
        # TODO: if there are no successors (b/c of questioned student) modify the if 13j and return case
        # Maybe: also add an extra variable in current node so we precisely say it's a "waiting" node generation
        successors = []
        initialPosition = currentNode.position
        for i in range(len(self.movesX)):
            newPosition = {
                    "x": initialPosition["x"] + self.movesX[i],
                    "y": initialPosition["y"] + self.movesY[i]
                    }
            if self.validPosition(initialPosition, newPosition):
                id = self.matrix[newPosition["x"]][newPosition["y"]]
                if id == "liber":
                    continue
                if (not currentNode.isInPath(id) \
                        and (not id in self.upset[currentNode.id]) \
                        and (not self.isInPerimeter(currentNode, newPosition))):
                    nextIndex = currentNode.index + 1
                    nextQuestioned = None
                    if nextIndex // 4 < len(self.questioned):
                        nextQuestioned = self.questioned[nextIndex // 4]
                    newNode = NodeTransition(
                            nextIndex,
                            id, 
                            nextQuestioned,
                            newPosition,
                            currentNode,
                            currentNode.g + 1
                            )
                    successors.append(newNode)
        return successors
    
    def isInPerimeter(self, currentNode, newPosition):
        """
        Checks if the new possible student is in the perimeter 
        of the currently questioned student
        returns the boolean answer
        """
        if currentNode.currentQuestioned != None:
            pos = self.getIndexNode(currentNode.currentQuestioned)
            if (pos["x"] - 1) <= newPosition["x"] and newPosition["x"] <= (pos["x"] + 1):
                yIndex = newPosition["y"] // 2
                if (pos["y"] // 2 - 1) <= yIndex and yIndex <= (pos["y"] // 2 + 1):
                    return True
        return False


    def validPosition(self,pastPosition, newPosition):
        """
        Checks if the new position is a valid move (inside the matrix)
        and checks if it moves to another 2k column, to be in the last 2 rows
        returns the boolean answer
        """
        if newPosition["x"] >= len(self.matrix) or newPosition["x"] < 0:
            return False
        if newPosition["y"] >= len(self.matrix[newPosition["x"]]) or newPosition["y"] < 0:
            return False
        if newPosition["y"] // 2 != pastPosition["y"] // 2 and newPosition["x"] < (len(self.matrix) - 2):
            return False
        return True

def uniform_cost_search(graph):
    que = [NodeTransition(0, graph.start, None, graph.getIndexNode(graph.start), None, 0)]

    while len(que) > 0:
        print("Current queue:")
        print(que)
        currentNode = que.pop(0)

        if graph.test_goal(currentNode):
            print("---------------------------------------------------")
            print("Solution: ", end = "")
            print(currentNode.showPath())
            print ("Cost: {}".format(currentNode.g))
            return

        successors = graph.generateSuccessors(currentNode)
        for node in successors:
            i = 0
            found_order = False
            while i < len(que):
                if que[i].g > node.g:
                    found_order = True
                    break
                i += 1
            if found_order:
                que.insert(i, node)
            else:
                que.append(node)
            


def main():
    data = read_file("./input/1_in.in")
    # we still need to do some data precomputations
    upset = {}
    nodes = {}
    for i in range(len(data["data"])):
        for j in range(len(data["data"][i])):
            if data["data"][i][j] != "liber":
                upset[data["data"][i][j]] = []
                nodes[data["data"][i][j]] = {"x": i, "y": j}
    for elem in data["upset"]:
        upset[elem[0]].append(elem[1])
        upset[elem[1]].append(elem[0])
    timeQuestioned = data["questioned"][0]
    questioned = list(data["questioned"][1:])

    graph = Graph(
            copy.deepcopy(nodes),
            copy.deepcopy(data["data"]),
            data["message"][0],
            data["message"][1],
            questioned,
            upset,
            timeQuestioned
            )

    uniform_cost_search(graph)

if __name__ == "__main__":
    main()

