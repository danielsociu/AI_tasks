import time
import random
import re
import math
import sys
import os
import copy

flatten = lambda t: [item for sublist in t for item in sublist]
printGlobally = False

def read_file(path):
    """
    Reads file at path, the output will be a dictionary with:
    {   data: the matrix of data,
        upset: the list of pairs of upset kids
        questioned: the list of quoestioned kids
        message: the starting and goal point of message}
    Returns:
        The data from the files
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


class NodeTransition:
    def __init__(self, index, id, currentQuestioned, position, father, cost, h, around):
        self.index = index
        self.id = id
        self.currentQuestioned = currentQuestioned
        self.position = position
        self.father = father
        self.g = cost
        self.h = h
        self.f = self.g + self.h
        self.around = around

    def getPath(self):
        path = [self.id]
        node = self.father
        while node is not None:
            path.insert(0, node.id)
            node = node.father
        return path

    def showPath(self):
        strPath = str(self.index + 1) + "." + self.id
        node = self.father
        prev = self
        while node is not None:
            if prev.position["y"] < node.position["y"]:
                if prev.position["y"] // 2 == node.position["y"] // 2:
                    strPath =  " < " + strPath
                elif prev.position["y"] // 2 < node.position["y"] // 2:
                    strPath =  " << " + strPath
            elif prev.position["y"] > node.position["y"]:
                if prev.position["y"] // 2 == node.position["y"] // 2:
                    strPath =  " > " + strPath
                elif prev.position["y"] // 2 > node.position["y"] // 2:
                    strPath =  " >> " + strPath
            elif prev.position["x"] < node.position["x"]:
                strPath =  " ^ " + strPath
            else: 
                strPath =  " v " + strPath
            strPath = str(node.index + 1) + "." + node.id + strPath
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
        if self.currentQuestioned != None:
            myStr += self.currentQuestioned
        myStr += ", cost = " + str(self.g) + "), |||| "
        return myStr



class Graph:
    movesX = [0, 0, 1, -1]
    movesY = [1, -1, 0, 0]
    def __init__(self, nodes, matrix, start, goal, questioned, upset, timeQuestioned, nameHeuristic = "naive"):
        self.nodes = copy.deepcopy(nodes)
        self.matrix = copy.deepcopy(matrix)
        self.start = start
        self.goal = goal
        self.questioned = copy.deepcopy(questioned)
        self.upset = copy.deepcopy(upset)
        self.timeQuestioned = timeQuestioned
        self.nameHeuristic = nameHeuristic

    def getIndexNode(self, name):
        """
        self.nodes will have the format node: location {x: , y: }
        since every node-name is unique
        """
        return self.nodes[name]

    def check_upset(self, id1, id2):
        if id1 == 'liber' or id2 == 'liber':
            return True
        if id2 in self.upset[id1]:
            return True
        return False

    def check_solvable(self):
        posStart = self.getIndexNode(self.start)
        posGoal = self.getIndexNode(self.goal)
        starty = (min(posStart["y"], posGoal['y']) // 2) * 2 
        endy = (max(posStart["y"], posGoal['y']) // 2) * 2
        for i in range(starty + 1, endy, 2):
            length = len(self.matrix) - 1
            if self.check_upset(self.matrix[length][i], self.matrix[length][i + 1])\
                    and self.check_upset(self.matrix[length - 1][i], self.matrix[length - 1][i + 1]):
                return False
        if starty != endy:
            for i in range(posStart['x'], len(self.matrix) - 2):
                j = posStart['y'] - (posStart['y'] % 2)
                if self.check_upset(self.matrix[i][j], self.matrix[i + 1][j]) \
                        and self.check_upset(self.matrix[i][j + 1], self.matrix[i + 1][j + 1]):
                    return False
            for i in range(posGoal['x'], len(self.matrix) - 2):
                j = posGoal['y'] - (posGoal['y'] % 2)
                if self.check_upset(self.matrix[i][j], self.matrix[i + 1][j]) \
                        and self.check_upset(self.matrix[i][j + 1], self.matrix[i + 1][j + 1]):
                    return False
        return True

    def test_goal(self, currentNode):
        return self.goal == currentNode.id

    def generateSuccessors(self, currentNode):
        # TODO: if there are no successors (b/c of questioned student) modify the if 13j and return case
        # Maybe: also add an extra variable in current node so we precisely say it's a "waiting" node generation
        successors = []
        successors_busy = []
        initialPosition = currentNode.position
        inPerimeterOfQuestioned = 0
        for i in range(len(self.movesX)):
            newPosition = {
                    "x": initialPosition["x"] + self.movesX[i],
                    "y": initialPosition["y"] + self.movesY[i]
                    }
            if self.validPosition(initialPosition, newPosition):
                id = self.matrix[newPosition["x"]][newPosition["y"]]
                if id == "liber":
                    continue
                if self.isInPerimeter(currentNode, newPosition):
                    inPerimeterOfQuestioned += 1
                    continue
                if (not id in self.upset[currentNode.id]) \
                        and (not self.isInPerimeter(currentNode, newPosition)):
                    nextIndex = currentNode.index + 1

                    nextQuestioned = None
                    if nextIndex // self.timeQuestioned < len(self.questioned):
                        nextQuestioned = self.questioned[nextIndex // self.timeQuestioned]

                    newCost = currentNode.g
                    if self.movesX[i] != 0:
                        newCost += 1
                    elif newPosition["y"] // 2 != initialPosition["y"] // 2:
                        newCost += 2

                    if (not currentNode.isInPath(id)): 
                        newNode = NodeTransition(
                                nextIndex,
                                id, 
                                nextQuestioned,
                                newPosition,
                                currentNode,
                                newCost,
                                self.calculate_h(id),
                                False
                                )
                        successors.append(newNode)
                    elif currentNode.currentQuestioned != None:
                        newNode = NodeTransition(
                                nextIndex,
                                id, 
                                nextQuestioned,
                                newPosition,
                                currentNode,
                                newCost,
                                self.calculate_h(id),
                                True
                                )
                        successors_busy.append(newNode) 
        # if we have a succesor that's been denied and no succesors we generate some other variants
        if len(successors) == 0 and inPerimeterOfQuestioned > 0:
            return successors_busy
        return successors

    def calculate_h(self, name):
        if self.nameHeuristic == "naive":
            return 1
        if self.goal == name:
            return 0
        if self.nameHeuristic == "admissible1":
            val = 0
            pos = self.getIndexNode(name)
            posGoal = self.getIndexNode(name)
            if posGoal["y"] // 2 != pos["y"] // 2:
                val += math.sqrt((posGoal["y"] - pos["y"]) ** 2 + (posGoal["x"] - pos["x"]) ** 2)
            else:
                val += abs(posGoal['x'] - pos['x'])
            return val
        if self.nameHeuristic == "admissible2":
            val = 0
            pos = self.getIndexNode(name)
            posGoal = self.getIndexNode(name)
            val += abs(posGoal["y"] // 2 - pos["y"] // 2) * 2
            if posGoal["y"] // 2 != pos["y"] // 2:
                val += 2*abs(posGoal['x'] - pos['x'])
            return val
        if self.nameHeuristic == "inadmissible":
            val = 0
            pos = self.getIndexNode(name)
            val += pos['x'] * 2
            val += pos['y']
            return val


    
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
                if ((pos["y"] + 1) // 2) <= yIndex and yIndex <= ((pos["y"] - 1) // 2):
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

def uniform_cost_search(graph, numberOfSolutions, fout, startTime, timeoutTime):
    timeout = time.time()
    total_generated = 1
    max_generated = 1
    firstQuestioned = None
    if graph.check_solvable() == False:
        fout.write("~~~~~~~~~~~~~~ No solutions ~~~~~~~~~~~~~~~\n")
    if len(graph.questioned) > 0: 
        firstQuestioned = graph.questioned[0]
    que = [
            NodeTransition(
                0, 
                graph.start,
                firstQuestioned, 
                graph.getIndexNode(graph.start), 
                None, 
                0, 
                0, 
                False
            )]

    while len(que) > 0:
        if printGlobally:
            print("Current queue:")
            print(que)
        currentNode = que.pop(0)

        if graph.test_goal(currentNode):
            write_to_file(fout, currentNode, startTime, max_generated, total_generated)
            max_generated = len(que)
            if printGlobally:
                print("---------------------------------------------------")
                print("Solution: ", end = "")
                print(currentNode.showPath())
                print ("Cost: {}".format(currentNode.g))
            numberOfSolutions -= 1
            # input()
            if numberOfSolutions == 0:
                return

        if timeoutTime < (round(1000*(time.time() - timeout))):
            if printGlobally:
                print("Timeout exceeded")
            fout.write("~~~~~~~~~~~~~~ Timeout exceeded ~~~~~~~~~~~~~~~\n")
            return

        successors = graph.generateSuccessors(currentNode)
        max_generated = max(max_generated, len(successors) + len(que))
        total_generated += len(successors)
        for node in successors:
            i = 0
            found_order = False
            while i < len(que):
                if que[i].g >= node.g:
                    found_order = True
                    break
                i += 1
            if found_order:
                que.insert(i, node)
            else:
                que.append(node)

def a_star(graph, numberOfSolutions, fout, startTime, timeoutTime):
    timeout = time.time()
    total_generated = 1
    max_generated = 1
    firstQuestioned = None
    if graph.check_solvable() == False:
        fout.write("~~~~~~~~~~~~~~ No solutions ~~~~~~~~~~~~~~~\n")
    if len(graph.questioned) > 0: 
        firstQuestioned = graph.questioned[0]
    que = [
            NodeTransition(
                0, 
                graph.start,
                firstQuestioned, 
                graph.getIndexNode(graph.start), 
                None, 
                0, 
                graph.calculate_h(graph.start), 
                False
            )]

    while len(que) > 0:
        if printGlobally: 
            print("Current queue:")
            print(que)
        currentNode = que.pop(0)

        if graph.test_goal(currentNode):
            write_to_file(fout, currentNode, startTime, max_generated, total_generated)
            max_generated = len(que)
            if printGlobally: 
                print("---------------------------------------------------")
                print("Solution: ", end = "")
                print(currentNode.showPath())
                print ("Cost: {}".format(currentNode.g))
            numberOfSolutions -= 1
            # input()
            if numberOfSolutions == 0:
                return

        if timeoutTime < (round(1000*(time.time() - timeout))):
            if printGlobally: 
                print("Timeout exceeded")
            fout.write("~~~~~~~~~~~~~~ Timeout exceeded ~~~~~~~~~~~~~~~\n")
            return

        successors = graph.generateSuccessors(currentNode)
        max_generated = max(max_generated, len(successors) + len(que))
        total_generated += len(successors)
        for node in successors:
            i = 0
            found_order = False
            while i < len(que):
                if que[i].f > node.f or (que[i].f == node.f and que[i].g <= node.g):
                    found_order = True
                    break
                i += 1
            if found_order:
                que.insert(i, node)
            else:
                que.append(node)

def a_star_optimised(graph, fout, startTime, timeoutTime):
    timeout = time.time()
    total_generated = 1
    max_generated = 1
    firstQuestioned = None
    if graph.check_solvable() == False:
        fout.write("~~~~~~~~~~~~~~ No solutions ~~~~~~~~~~~~~~~\n")
    if len(graph.questioned) > 0: 
        firstQuestioned = graph.questioned[0]
    que = [
            NodeTransition(
                0, 
                graph.start,
                firstQuestioned, 
                graph.getIndexNode(graph.start), 
                None, 
                0, 
                graph.calculate_h(graph.start),
                False
            )]

    que_closed = []
    while len(que) > 0:
        if printGlobally: 
            print("Current queue:")
            print(que)
        currentNode = que.pop(0)
        que_closed.append(currentNode)

        if graph.test_goal(currentNode):
            write_to_file(fout, currentNode, startTime, max_generated, total_generated)
            max_generated = len(que)
            if printGlobally: 
                print("---------------------------------------------------")
                print("Solution: ", end = "")
                print(currentNode.showPath())
                print ("Cost: {}".format(currentNode.g))
            # input()
            return

        if timeoutTime < (round(1000*(time.time() - timeout))):
            if printGlobally: 
                print("Timeout exceeded")
            fout.write("~~~~~~~~~~~~~~ Timeout exceeded ~~~~~~~~~~~~~~~\n")
            return

        successors = graph.generateSuccessors(currentNode)
        
        max_generated = max(max_generated, len(successors) + len(que))
        total_generated += len(successors)
        # Optimisations
        for s in successors:
            found = False
            for node in que:
                if s.id == node.id:
                    found = True
                    if s.f > node.f and s.around == False:
                        successors.remove(s)
                    elif node.around == False:
                        que.remove(node)
                    break
            if not found:
                for node in que_closed:
                    if s.id == node.id:
                        if s.f > node.f and s.around == False:
                            successors.remove(s)
                        elif node.around == False:
                            que_closed.remove(node)
                        break

        for node in successors:
            i = 0
            found_order = False
            while i < len(que):
                if que[i].f > node.f or (que[i].f == node.f and que[i].g <= node.g):
                    found_order = True
                    break
                i += 1
            if found_order:
                que.insert(i, node)
            else:
                que.append(node)

def IDA_star(graph, numberOfSolutions, fout, startTime, timeoutTime):
    timeout = time.time()
    total_generated = 1
    firstQuestioned = None
    if graph.check_solvable() == False:
        fout.write("~~~~~~~~~~~~~~ No solutions ~~~~~~~~~~~~~~~\n")
    if len(graph.questioned) > 0: 
        firstQuestioned = graph.questioned[0]
    initialNode = NodeTransition(
                0, 
                graph.start,
                firstQuestioned, 
                graph.getIndexNode(graph.start), 
                None, 
                0, 
                graph.calculate_h(graph.start),
                False
            )
    limit = initialNode.f

    while True:
        if printGlobally: 
            print("Limit of tree depth: ", limit)
        numberOfSolutions, result, generated = build_tree(graph, initialNode, limit, numberOfSolutions, fout, startTime, timeout, timeoutTime, 0, total_generated)
        total_generated += generated
        if result == "done":
            break
        if result == float("inf"):
            if printGlobally: 
                print("No solutions found")
            break

        if timeoutTime < (round(1000*(time.time() - timeout))) or result == "timeout":
            if printGlobally: 
                print("Timeout exceeded")
            fout.write("~~~~~~~~~~~~~~ Timeout exceeded ~~~~~~~~~~~~~~~\n")
            return

        limit = result
        if printGlobally: 
            print(">>>>> New limit: {} <<<<<<".format(limit))
        # input()

def build_tree(graph, currentNode, limit, numberOfSolutions, fout, startTime, timeout, timeoutTime,currently_generated, total_generated):
    if printGlobally: 
        print("Reached: ", currentNode)
        print("G and H and F: {} {} {}".format(currentNode.g, currentNode.h, currentNode.f))
        print("test_goal: ", graph.test_goal(currentNode))
    currentlyGen = 0
    if currentNode.f > limit:
        return numberOfSolutions, currentNode.f, 0
    if graph.test_goal(currentNode) and currentNode.f == limit:
        write_to_file(fout, currentNode, startTime, currently_generated, total_generated)
        if printGlobally: 
            print("---------------------------------------------------")
            print("Solution: ", end = "")
            print(currentNode.showPath())
            print("Limit: {}".format(limit))
            print ("Cost: {}".format(currentNode.g))
        # input()
        numberOfSolutions -= 1
        if numberOfSolutions == 0:
            return 0, "done", 0
    successors = graph.generateSuccessors(currentNode)
    minim = float("inf")

    if timeoutTime < (round(1000*(time.time() - timeout))):
        return 0, "timeout", 0


    for s in successors:
        numberOfSolutions, result, generated = build_tree(graph, s, limit, numberOfSolutions, fout, startTime, timeout, timeoutTime, max(currently_generated, currentlyGen), total_generated + currentlyGen)
        currentlyGen += generated
        if result == "done":
            return 0, "done", currentlyGen
        if result == "timeout":
            return 0, "timeout", currentlyGen
        if printGlobally: 
            print("Comparing ", result, " with ", minim)
        if result < minim:
            minim = result
            if printGlobally: 
                print("New minim: ", minim)
    currentlyGen += len(successors)
    return numberOfSolutions, minim, currentlyGen

def write_to_file(fout, currentNode, startTime, maxGenerated, totalGenerated):
    fout.write("Solution: ")
    fout.write(currentNode.showPath() + '\n')
    fout.write("Length: {}\n".format(currentNode.index + 1))
    fout.write("Cost: {}\n".format(currentNode.g))
    fout.write("Time: {}\n".format(round(1000*(time.time() - startTime))))
    fout.write("Max nodes: {}\n".format(maxGenerated))
    fout.write("Total nodes: {}\n".format(totalGenerated))
    fout.write("---------------------------------------------------\n")


def prelucrate_data(inData):
    data = copy.deepcopy(inData)
    upset = {}
    nodes = {}
    for i in range(len(data["data"])):
        for j in range(len(data["data"][i])):
            if data["data"][i][j] != "liber":
                upset[data["data"][i][j]] = []
                nodes[data["data"][i][j]] = {"x": i, "y": j}
    for elem in data["upset"]:
        if elem[0] not in upset.keys():
            upset[elem[0]] = [elem[1]]
        else:
            upset[elem[0]].append(elem[1])
        if elem[1] not in upset.keys():
            upset[elem[1]] = [elem[0]]
        else:
            upset[elem[1]].append(elem[0])
    timeQuestioned = int(data["questioned"][0])
    questioned = list(data["questioned"][1:])
    return nodes, data["data"], data["message"][0], data["message"][1], questioned, upset, timeQuestioned

def main():
    if (len(sys.argv) < 4):
        print("Wrong number of arguments!")
        return
    folderIn = sys.argv[1]
    folderOut = sys.argv[2]
    NSol = int(sys.argv[3])
    timeoutTime = int(sys.argv[4])

    startTimer = time.time()

    # Preparing data out
    if not os.path.exists(folderOut):
        os.mkdir(folderOut)

    # Data in
    listFolder = os.listdir(folderIn)
    for file in listFolder:
        fileName = file.split(".")[0]
        foutUCS = open("./{}/{}_1_UCS.out".format(folderOut, fileName), "w")
        data = read_file("./{}/{}".format(folderIn, file))
        nodes, matrix, start, goal, questioned, upset, timeQuestioned = prelucrate_data(data)
        graphUCS = Graph(
                nodes,
                matrix,
                start,
                goal,
                questioned,
                upset,
                timeQuestioned,
                "naive"
                )
        foutUCS.write(">>>>>> UCS solutions <<<<<<<<\n")
        uniform_cost_search(graphUCS, NSol, foutUCS, startTimer, timeoutTime)

        heuristics = ["naive", "admissible1", "admissible2", "inadmissible"]
        fileWriting = "w"

        for heurisitic in heuristics:
            graphAS = Graph(
                    nodes,
                    matrix,
                    start,
                    goal,
                    questioned,
                    upset,
                    timeQuestioned,
                    heurisitic
                    )
            graphASO = Graph(
                    nodes,
                    matrix,
                    start,
                    goal,
                    questioned,
                    upset,
                    timeQuestioned,
                    heurisitic
                    )
            graphIDA = Graph(
                    nodes,
                    matrix,
                    start,
                    goal,
                    questioned,
                    upset,
                    timeQuestioned,
                    heurisitic
                    )
            foutAS = open("./{}/{}_2_AS.out".format(folderOut, fileName), fileWriting)
            foutAS.write("\n>>>>>> A star solutions: " + heurisitic + " <<<<<<<<\n")
            a_star(graphAS, NSol, foutAS, startTimer, timeoutTime)
            foutASO = open("./{}/{}_3_ASO.out".format(folderOut, fileName), fileWriting)
            foutASO.write("\n>>>>>> A star optimised solutions: " + heurisitic + " <<<<<<<<\n")
            a_star_optimised(graphASO, foutASO, startTimer, timeoutTime)
            foutIDA = open("./{}/{}_4_IDA.out".format(folderOut, fileName), fileWriting)
            foutIDA.write("\n>>>>>> IDA solutions: " + heurisitic + " <<<<<<<<\n")
            IDA_star(graphIDA, NSol, foutIDA, startTimer, timeoutTime)
            fileWriting = "a"
    # we still need to do some data precomputations

    
if __name__ == "__main__":
    main()

