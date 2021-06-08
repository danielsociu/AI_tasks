import time
import copy
import json
import pygame
import sys

class Cell:
    """
    A cell of the matrix
    """
    wallWidth = 15
    backgroundCell = (255,255,255)
    lineColor = (33, 222, 71)
    display = None

    def __init__(self, left = None, top = None, width = None, height = None, line = None, column = None, code = 0, wall = None):
        if width != None and height != None:
            self.rectangle = pygame.Rect(left, top, width, height)
        self.wall = [None, None, None, None]
        self.code = code
        self.line = line
        self.width = width
        self.height = height
        self.column = column
        self.top = top
        self.left = left
        if wall:
            if wall != "skip":
                for i in len(wall):
                    self.wall[i] = wall[i]
        else:
            if line > 0:
                self.wall[0] = pygame.Rect(left, top - 1 - self.__class__.wallWidth // 2,\
                        width, self.__class__.wallWidth)
            else:
                self.code += 2**0
            if column < Interface.columnsNumber - 1:
                self.wall[1] = pygame.Rect(left + width - self.__class__.wallWidth // 2, top,\
                        self.__class__.wallWidth, height)
            else:
                self.code += 2**1
            if line < Interface.linesNumber - 1:
                self.wall[2] = pygame.Rect(left, top + height - self.__class__.wallWidth // 2, \
                        width, self.__class__.wallWidth)
            else:
                self.code += 2**2
            if column > 0:
                self.wall[3] = pygame.Rect(left - self.__class__.wallWidth // 2, top, \
                        self.__class__.wallWidth, height)
            else:
                self.code += 2**3

    def copyCell(self):
        """
        copy a cell for our min_max and alpha_beta algorithms
        copies only what's necessary for finding best solution
        """
        newCell = Cell(
                column = self.column,
                line = self.line,
                code = self.code,
                wall = "skip"
                )
        return newCell

    def drawCell(self, index = None, background = None):
        """
        Draws a cell to the display
        """
        if background:
            pygame.draw.rect(Cell.display, background, self.rectangle)
        else:
            pygame.draw.rect(Cell.display, self.__class__.backgroundCell, self.rectangle)
        byte = 1
        for i in range(4):
            if self.code & byte or index == i:
                if self.wall[i]:
                    pygame.draw.rect(Cell.display, self.__class__.lineColor, self.wall[i])
            byte *= 2

class Piece:
    """
    A game piece (red/blue)
    """
    def __init__(self, image, moves, name, dimension, line, column, goal):
        self.image = image 
        self.name = name
        self.dimension = dimension
        self.line = line
        self.column = column
        self.goal = goal
        self.moves = moves
    def copyPiece(self):
        newPiece = Piece(None, self.moves, self.name, self.dimension, self.line, self.column, self.goal)
        return newPiece
    def __repr__(self):
        string = self.name + ' at: ({},{}) moves left: {}'.format(self.line, self.column, self.moves)
        return string


class Interface:
    """
    Game interface, it has data about drawing and
    the cells of the game (walls as well)
    """
    screenColor = (0, 0, 0)
    screen = None
    linesNumber = None
    columnsNumber = None 
    cellWidth = None
    cellHeight = None
    cellPadding = None
    imageDimension = None

    def __init__(self, cellMatrix = None):
        if cellMatrix:
            self.cellMatrix = cellMatrix
        else:
            self.cellMatrix = [[
                Cell(
                    left = col * (self.cellWidth + 1),
                    top = line * (self.cellHeight + 1),
                    width = self.cellWidth,
                    height = self.cellHeight,
                    line = line,
                    column = col
                    ) for col in range(Interface.columnsNumber)]
                for line in range(Interface.linesNumber)
                ]

    def copyMatrix(self):
        """
        copy a matrix with minimal needs
        """
        cellMatrix = [[
                cell.copyCell() for cell in line]
            for line in self.cellMatrix
            ]
        return cellMatrix

    @classmethod
    def initialize(
            cls,
            linesNumber = 9,
            columnsNumber = 9,
            cellWidth = 100,
            cellHeight = 100,
            cellPadding = 8,
            screenColor = (0, 0, 0),
            screen = None
            ):
        """
        Initializes the class defaults
        """
        cls.linesNumber = linesNumber
        cls.columnsNumber = columnsNumber
        cls.cellWidth = cellWidth
        cls.cellHeight = cellHeight
        cls.cellPadding = cellPadding
        cls.imageDimension = min(cellWidth, cellHeight) - 2 * cellPadding
        cls.screen = screen

    def drawImage(self, image, cell):
        """
        Draws an image to screen
        """
        self.screen.blit(
                image,
                (
                    cell.rectangle.left + self.cellPadding,
                    cell.rectangle.top + self.cellPadding
                ))

    #(cell, k, wall, i, j)
    """
    Draws all the screen, with cells and pieces
    """
    def drawGameScreen(self, pieces, walls = None, highlightCells = None):
        Interface.screen.fill(self.screenColor)
        for i, line in enumerate(self.cellMatrix):
            for j, cell in enumerate(line):
                if highlightCells and cell in highlightCells:
                    cell.drawCell(background = (255,255,0))
                else:
                    cell.drawCell()

        for piece in pieces.values():
            self.drawImage(piece.image, self.cellMatrix[piece.line][piece.column])
        if walls:
            for wall in walls:
                wall[0].drawCell(index = wall[1])
        pygame.display.update()

class Game:
    GMIN = None
    GMAX = None
    max_score = 0

    def __init__(self, interface, red, blue):
        self.red = red
        self.blue = blue
        self.pieces = {"red": self.red, "blue": self.blue}
        self.interface = interface
        self.__class__.max_score = Interface.linesNumber * Interface.columnsNumber
    
    def drawGameScreen(self, **kwargs):
        self.interface.drawGameScreen(self.pieces, **kwargs)

    def final(self):
        """
        Check if the game is in a final state
        """
        if self.red.line == self.red.goal:
            return "red"
        elif self.blue.line == self.blue.goal:
            return "blue"
        return False
    
    def estimateScore(self, depth):
        """
        Estimation score, it considers the number of moves to 
        to goal and how many walls the user has left
        """
        final = self.final()
        if final == self.__class__.GMAX:
            return self.__class__.max_score + depth
        elif final == self.__class__.GMIN:
            return -(self.__class__.max_score + depth)
        else:
            gmaxPlayer = self.pieces[self.__class__.GMAX]
            gminPlayer = self.pieces[self.__class__.GMIN]
            # print ((self.semiLee(self.__class__.GMIN) , self.semiLee(self.__class__.GMAX)))
            return (self.semiLee(self.__class__.GMIN) - self.semiLee(self.__class__.GMAX) + 
                    (gmaxPlayer.moves - gminPlayer.moves)
                )

    def estimateScore2(self, depth):
        """
        Estimation score that checks how many lines till solution
        """
        final = self.final()
        if final == self.__class__.GMAX:
            return self.__class__.max_score + depth
        elif final == self.__class__.GMIN:
            return -(self.__class__.max_score + depth)
        else:
            gmaxPlayer = self.pieces[self.__class__.GMAX]
            gminPlayer = self.pieces[self.__class__.GMIN]
            return (abs(gminPlayer.line - gminPlayer.goal) - abs(gmaxPlayer.line - gmaxPlayer.goal) +
                    (gmaxPlayer.moves - gminPlayer.moves)
                    )

    @classmethod
    def getOppositePlayer(cls, player):
        return cls.GMAX if player == cls.GMIN else cls.GMIN

    def getOppositeWalls(self, x, y, index):
        """
        Gets the opposite wall to verify a cut,
        """
        nextIndex = None
        if x == 1:
            nextIndex = (index + (1 if (index == 1) else -1)) % 4 
        elif x == -1:
            nextIndex = (index + (1 if (index == 3) else -1)) % 4 
        elif y == 1:
            nextIndex = (index + (1 if (index == 0) else -1)) % 4 
        else:
            nextIndex = (index + (1 if (index == 2) else -1)) % 4 
        return nextIndex

    def checkOppositeWall(self,cell1, cell2, x, y, index):
        """
        checks if the wall we want to build will cut another one
        """
        nextIndex1 = self.getOppositeWalls(x, y, index)
        nextIndex2 = self.getOppositeWalls(x, y, (index + 2) % 4)
        if cell1.code & 2**nextIndex1 and cell2.code & 2**nextIndex2:
            return False
        return True


    def checkValidationWall(self, i, j, moveX, moveY, index):
        """
        Checks if the wall is valid (if there no other wall there)
        """
        nextCell = self.interface.cellMatrix[i + moveX][j + moveY]
        if (nextCell.code & 2**index):
            return None
        return (nextCell, index, nextCell.wall[index], i + moveX, j + moveY)

    def getWallContinuation(self, currentWall):
        """
        Gets the wall continuation, when we send 2 walls
        Priority will be towards + 1 in either direction
        """
        affectedCells = []
        bad = 0
        if (currentWall[0][1] & 1) != (currentWall[1][1] & 1) or (currentWall[0][0].code & 2**currentWall[0][1]):
            return currentWall
        for other, (cell, index, wall, i, j) in enumerate(currentWall):
            nextTuple = None
            otherCell = currentWall[other ^ 1][0]
            if index & 1:
                if i < Interface.linesNumber - 1 and self.checkOppositeWall(cell, otherCell, 1, 0, index):
                    move = 1
                    otherTuple = currentWall[other ^ 1]
                    nextTuple = self.checkValidationWall(i, j, move, 0, index)
                if not nextTuple and i > 0 and self.checkOppositeWall(cell, otherCell, -1, 0, index):
                    move = -1
                    nextTuple = self.checkValidationWall(i, j, move, 0, index)
            else:
                if j < Interface.columnsNumber - 1 and self.checkOppositeWall(cell, otherCell, 0, 1, index):
                    move = 1
                    nextTuple = self.checkValidationWall(i, j, 0, move, index)
                if not nextTuple and j > 0 and self.checkOppositeWall(cell, otherCell, 0, -1, index):
                    move = -1
                    nextTuple = self.checkValidationWall(i, j, 0, move, index)

            if nextTuple != None:
                affectedCells.append(nextTuple)
        affectedCells += currentWall
        return affectedCells

    def insideCheck(self, pos):
        """
        Checks if the position dictionary is inside the matrix
        """
        if 0 <= pos['x'] < Interface.linesNumber and 0 <= pos['y'] < Interface.columnsNumber:
            return True
        return False

    def semiLee(self, playerName):
        """
        Lee that gets the shortest path to goal, if it exists
        or 0 otherwise
        """
        nextX = [-1, 0, 1, 0]
        nextY = [0, 1, 0, -1]
        player = self.pieces[playerName]
        visited = [[0] * Interface.columnsNumber for x in range(Interface.linesNumber)]
        # print(visited)
        que = []
        que.append({'x': player.line, 'y': player.column})
        visited[player.line][player.column] = 1
        goal = player.goal
        if player.line == goal:
            return 0
        while len(que) != 0:
            pos = que.pop(0)
            cell = self.interface.cellMatrix[pos['x']][pos['y']]
            for index in range(len(nextX)):
                nextPos = {
                    'x': pos['x'] + nextX[index],
                    'y': pos['y'] + nextY[index]
                }
                if not (cell.code & (2 ** index)) and self.insideCheck(nextPos) and not visited[nextPos['x']][nextPos['y']]:
                    if nextPos['x'] == goal:
                        return visited[pos['x']][pos['y']] 
                    que.append(nextPos)
                    visited[nextPos['x']][nextPos['y']] = visited[pos['x']][pos['y']] + 1
        return 0
    
    def checkWallBlock(self, affectedCells):
        """
        Checks if adding a wall would block any of the pieces
        and prevent them from getting to their goal
        """
        ok = True
        for (cell, index, wall, i, j) in affectedCells:
            cell.code |= 2 ** index
        if self.semiLee(self.blue.name) and self.semiLee(self.red.name):
            ok = False
        # print (ok)
        for (cell, index, wall, i, j) in affectedCells:
            cell.code -= 2 ** index
        return ok

    def getNextCells(self, playerName):
        """
        Gets next cells that the player can move to
        Also checks if there is another player nearby or a wall block
        """
        nextX = [-1, 0, 1, 0]
        nextY = [0, 1, 0, -1]
        player = self.pieces[playerName]
        currentCell = self.interface.cellMatrix[player.line][player.column]
        oppositePlayer = self.pieces[self.getOppositePlayer(playerName)]
        nextCells = []
        nextCells2 = []
        for index in range(len(nextX)):
            nextPos = {
                'x': player.line + nextX[index],
                'y': player.column + nextY[index]
            }
            if self.insideCheck(nextPos) and not(currentCell.code & (2 ** index)):
                curCell = self.interface.cellMatrix[nextPos['x']][nextPos['y']]
                if nextPos['x'] == oppositePlayer.line and nextPos['y'] == oppositePlayer.column:
                    oppositeCell = self.interface.cellMatrix[oppositePlayer.line][oppositePlayer.column]
                    for index2 in range(len(nextX)):
                        nextPos2 = {
                            'x': oppositePlayer.line + nextX[index2],
                            'y': oppositePlayer.column + nextY[index2]
                        }
                        if self.insideCheck(nextPos2) and not (oppositeCell.code & (2 ** index2))\
                                and not (nextPos2['x'] == player.line and nextPos2['y'] == player.column):
                            curCell2 = self.interface.cellMatrix[nextPos2['x']][nextPos2['y']]
                            nextCells2.append(curCell2)
                            if index == index2:
                                nextCells2 = [curCell2]
                                break
                else:
                    nextCells.append(curCell)
        return nextCells + nextCells2
    def wallMoves(self, player):
        """
        Generates all the wall moves (takes too much time to run)
        """
        games = []
        playerData = self.pieces[player]
        for i, line in enumerate(self.interface.cellMatrix):
            for j, cell in enumerate(line):
                if i < (Interface.linesNumber - 1):
                    curWall = [(cell, 2, cell.wall[2], i, j)]
                    cell2 = self.interface.cellMatrix[i + 1][j]
                    curWall.append((cell2, 0, cell2.wall[0], i + 1, j))
                    getWall = self.getWallContinuation(curWall)
                    if len(getWall) == 4 and not self.checkWallBlock(getWall):
                        cellMatrix = self.interface.copyMatrix()
                        # cellMatrix = copy.deepcopy(self.interface.cellMatrix)
                        for (myCell, myIndex, myWall, myI, myJ) in getWall:
                            cellMatrix[myI][myJ].code |= 2**myIndex
                        copyPlayerData = playerData.copyPiece()
                        copyPlayerData.moves -= 1
                        if player == 'red':
                            games.append(Game(Interface(cellMatrix), copyPlayerData, self.blue.copyPiece()))
                        else:
                            games.append(Game(Interface(cellMatrix), self.red.copyPiece(), copyPlayerData))
                if j < (Interface.columnsNumber - 1):
                    curWall = [(cell, 1, cell.wall[1], i, j)]
                    cell2 = self.interface.cellMatrix[i][j + 1]
                    curWall.append((cell2, 3, cell2.wall[3], i, j + 1))
                    getWall = self.getWallContinuation(curWall)
                    if  len(getWall) == 4 and not self.checkWallBlock(getWall):
                        cellMatrix = self.interface.copyMatrix()
                        # cellMatrix = copy.deepcopy(self.interface.cellMatrix)
                        for (myCell, myIndex, myWall, myI, myJ) in getWall:
                            cellMatrix[myI][myJ].code |= 2**myIndex
                        copyPlayerData = playerData.copyPiece()
                        copyPlayerData.moves -= 1
                        if player == 'red':
                            games.append(Game(Interface(cellMatrix), copyPlayerData, self.blue.copyPiece()))
                        else:
                            games.append(Game(Interface(cellMatrix), self.red.copyPiece(), copyPlayerData))
        return games

    def probableWallMoves(self, player):
        """
        Gets the relevant wall moves, the ones that are a continuation
        to an existing wall
        """
        games = []
        playerData = self.pieces[player]
        for i in range(len(self.interface.cellMatrix) - 1):
            line = self.interface.cellMatrix[i]
            for j in range(len(line) - 1):
                cell = line[j]
                # getting horizontal walls that are vertical continuations of others
                if cell.code & 4: 
                    jIndex = None
                    if j > 1 and not line[j - 1].code & 4:
                        jIndex = j - 1
                    elif j < len(line) - 1 and not line[j + 1].code & 4:
                        jIndex = j + 1
                    if jIndex == None:
                        continue
                    curWall = [
                            (line[jIndex], 2, cell.wall[2], i, jIndex),
                            (self.interface.cellMatrix[i + 1][jIndex], 0, cell.wall[0], i + 1, jIndex),
                            ]
                    getWall = self.getWallContinuation(curWall)
                    if len(getWall) == 4 and not self.checkWallBlock(getWall):
                        cellMatrix = self.interface.copyMatrix()
                        for (myCell, myIndex, myWall, myI, myJ) in getWall:
                            cellMatrix[myI][myJ].code |= 2**myIndex
                        copyPlayerData = playerData.copyPiece()
                        copyPlayerData.moves -= 1
                        if player == 'red':
                            games.append(Game(Interface(cellMatrix), copyPlayerData, self.blue.copyPiece()))
                        else:
                            games.append(Game(Interface(cellMatrix), self.red.copyPiece(), copyPlayerData))
                # getting vertical walls
                if cell.code & 2: 
                    iIndex = None
                    if i > 1 and not self.interface.cellMatrix[i - 1][j].code & 2:
                        iIndex = i - 1
                    elif i < len(self.interface.cellMatrix) - 1 and not self.interface.cellMatrix[i + 1][j].code & 2:
                        iIndex = i + 1
                    if iIndex == None:
                        continue
                    curWall = [
                            (self.interface.cellMatrix[iIndex][j], 1, cell.wall[1], iIndex, j),
                            (self.interface.cellMatrix[iIndex][j + 1], 3, cell.wall[3], iIndex, j + 1),
                            ]
                    getWall = self.getWallContinuation(curWall)
                    if len(getWall) == 4 and not self.checkWallBlock(getWall):
                        cellMatrix = self.interface.copyMatrix()
                        for (myCell, myIndex, myWall, myI, myJ) in getWall:
                            cellMatrix[myI][myJ].code |= 2**myIndex
                        copyPlayerData = playerData.copyPiece()
                        copyPlayerData.moves -= 1
                        if player == 'red':
                            games.append(Game(Interface(cellMatrix), copyPlayerData, self.blue.copyPiece()))
                        else:
                            games.append(Game(Interface(cellMatrix), self.red.copyPiece(), copyPlayerData))
        return games

    def userWalls(self, player):
        """
        Gets the wall moves around the opposite player
        """
        games = []
        oppositePlayer = self.pieces[self.getOppositePlayer(player)]
        playerData = self.pieces[player]
        moveX = [-1, 0, 1, 0]
        moveY = [0, 1, 0, -1]
        for index in range(len(moveX)):
            pos = {
                'x': oppositePlayer.line + moveX[index],
                'y': oppositePlayer.column + moveY[index],
            }
            if self.insideCheck(pos):
                cell = self.interface.cellMatrix[oppositePlayer.line][oppositePlayer.column]
                curWall = [
                        (cell, index, cell.wall[index], oppositePlayer.line, oppositePlayer.column),
                        (self.interface.cellMatrix[pos['x']][pos['y']], (index + 2) % 4, cell.wall[(index + 2) % 4], pos['x'], pos['y']),
                        ]
                getWall = self.getWallContinuation(curWall)
                if len(getWall) == 4 and not self.checkWallBlock(getWall):
                    cellMatrix = self.interface.copyMatrix()
                    for (myCell, myIndex, myWall, myI, myJ) in getWall:
                        cellMatrix[myI][myJ].code |= 2**myIndex
                    copyPlayerData = playerData.copyPiece()
                    copyPlayerData.moves -= 1
                    if player == 'red':
                        games.append(Game(Interface(cellMatrix), copyPlayerData, self.blue.copyPiece()))
                    else:
                        games.append(Game(Interface(cellMatrix), self.red.copyPiece(), copyPlayerData))
        return games


    def userMoves(self, player):
        """
        Gets the available user moves
        """
        games = []
        user_moves = self.getNextCells(player)
        for cell in user_moves:
            cellMatrix = self.interface.copyMatrix()
            newRed = self.red.copyPiece()
            newBlue = self.blue.copyPiece()
            if player == "red":
                newRed.line = cell.line
                newRed.column = cell.column
            else:
                newBlue.line = cell.line
                newBlue.column = cell.column
            games.append(Game(Interface(cellMatrix), newRed, newBlue))
        return games


class State:
    """
    A node in the graph for min_max and alpha_beta
    """
    def __init__(self, game, currentPlayer, depth, father = None, score = None):
        self.game = game
        self.currentPlayer = currentPlayer
        self.depth = depth
        self.father = father
        self.score = score

        self.possibleMoves = []
        self.chosenState = None

    def moves(self):
        """
        Generates all the possible moves for the current user
        """
        states = []
        if self.game.pieces[self.currentPlayer].moves > 0:
            wall_moves = self.game.probableWallMoves(self.currentPlayer)
            for game in wall_moves:
                states.append(State(game, self.game.getOppositePlayer(self.currentPlayer), self.depth - 1, father = self))
            user_walls = self.game.userWalls(self.currentPlayer)
            for game in user_walls:
                states.append(State(game, self.game.getOppositePlayer(self.currentPlayer), self.depth - 1, father = self))
        user_moves = self.game.userMoves(self.currentPlayer)
        for game in user_moves:
            states.append(State(game, self.game.getOppositePlayer(self.currentPlayer), self.depth - 1, father = self))
        return states

def min_max(state, scoring = 1):
    """
    min max algorithm to get the best move for the computer
    """
    if state.depth == 0  or state.game.final():
        if scoring == 1:
            state.score = state.game.estimateScore(state.depth)
        else:
            state.score = state.game.estimateScore2(state.depth)
        return state
    state.possibleMoves = state.moves()
    score_moves = [min_max(move) for move in state.possibleMoves]

    if state.currentPlayer == Game.GMAX:
        state.chosenState = max(score_moves, key = lambda x: x.score)
    else:
        state.chosenState = min(score_moves, key = lambda x: x.score)
    state.score = state.chosenState.score
    return state

def alpha_beta(alpha, beta, state, scoring = 1):
    """
    alpha beta algorithm to get the best move for the computer
    """
    if state.depth == 0  or state.game.final():
        if scoring == 1:
            state.score = state.game.estimateScore(state.depth)
        else:
            state.score = state.game.estimateScore2(state.depth)
        return state

    if alpha > beta:
        return state

    state.possibleMoves = state.moves()

    if state.currentPlayer == Game.GMAX:
        currentScore = float('-inf')
        for move in state.possibleMoves:
            newState = alpha_beta(alpha, beta, move)

            if (currentScore < newState.score):
                state.chosenState = newState
                currentScore = newState.score
            if (alpha < newState.score):
                alpha = newState.score
                if alpha >= beta:
                    break
    else:
        currentScore = float('inf')
        for move in state.possibleMoves:
            newState = alpha_beta(alpha, beta, move)

            if (currentScore > newState.score):
                state.chosenState = newState
                currentScore = newState.score
            if (beta > newState.score):
                beta = newState.score
                if alpha >= beta:
                    break
    state.score = state.chosenState.score
    return state


def byte_to_power(number):
    if number == 8:
        return 3
    elif number == 4:
        return 2
    elif number == 2:
        return 1
    else:
        return 0

class Button:
    def __init__(self, display = None, left = 0, top = 0, width = 0, height = 0, \
            backgroundColor = (53,80,115), selectedBackgroundColor = (89,134,194), text = "", font = "arial",\
            fontSize = 16, textColor = (255,255,255), value = ""):
        self.display = display
        self.left = left
        self.top = top
        self.width = width
        self.heght = height
        self.backgroundColor = backgroundColor
        self.selectedBackgroundColor = selectedBackgroundColor
        self.text = text
        self.font = font
        self.fontSize = fontSize
        self.textColor = textColor
        self.value = value
        fontObject = pygame.font.SysFont(self.font, self.fontSize)
        self.renderedText = fontObject.render(self.text, True, self.textColor)
        self.rectangle = pygame.Rect(left, top, width, height)
        self.rectangleText = self.renderedText.get_rect(center = self.rectangle.center)
        self.selected = False

    def selectButton(self, selected):
        self.selected = selected
        self.drawButton()

    def selectAfterCoord(self, coord):
        if self.rectangle.collidepoint(coord):
            self.selectButton(True)
            return True
        return False

    def updateRectangle(self):
        self.rectangle.left = self.left
        self.rectangle.top = self.top
        self.rectangleText = self.renderedText.get_rect(center = self.rectangle.center)

    def drawButton(self):
        color = self.selectedBackgroundColor if self.selected else self.backgroundColor
        pygame.draw.rect(self.display, color, self.rectangle)
        self.display.blit(self.renderedText, self.rectangleText)

class ButtonGroup:
    def __init__(self, buttonList = [], selectedIndex = 0, marginButtons = 10, left = 0, top = 0):
        self.buttonList = buttonList
        self.selectedIndex = selectedIndex
        self.buttonList[self.selectedIndex].selected = True
        self.marginButtons = marginButtons
        self.left = left
        self.top = top
        for button in self.buttonList:
            button.top = self.top
            button.left = left
            button.updateRectangle()
            left += (marginButtons + button.width)

    def selectAfterCoord(self, coord):
        for index, button in enumerate(self.buttonList):
            if button.selectAfterCoord(coord):
                self.buttonList[self.selectedIndex].selectButton(False)
                self.selectedIndex = index
                return True
        return False
    
    def drawButtons(self):
        for button in self.buttonList:
            button.drawButton()

    def getValue(self):
        return self.buttonList[self.selectedIndex].value

def draw_menu(display, game, size):
    buttons_algorithm = ButtonGroup(
            top = 140,
            left = size[1]//2 - 90,
            buttonList = [
                Button(display = display, width = 80, height = 30, text = "minimax", value = "minimax"),
                Button(display = display, width = 80, height = 30, text = "alphabeta", value = "alphabeta")
                ]
            )
    difficulty = ButtonGroup(
            top = 200,
            left = size[1]//2 - (80 * 3) / 2,
            buttonList = [
                Button(display = display, width = 80, height = 30, text = "beginner", value = 2),
                Button(display = display, width = 80, height = 30, text = "normal", value = 3),
                Button(display = display, width = 80, height = 30, text = "advanced", value = 4)
                ],
            selectedIndex = 1
            )
    buttons_player = ButtonGroup(
            top = 260,
            left = size[1]//2 - 60,
            buttonList = [
                Button(display = display, width = 50, height = 30, text = "red", value = "red"),
                Button(display = display, width = 50, height = 30, text = "blue", value = "blue")
                ]
            )
    play_mode = ButtonGroup(
            top = 320,
            left = size[1]//2 - (190 * 3) / 2,
            buttonList = [
                Button(display = display, width = 160, height = 30, text = "player vs player", value = "pp"),
                Button(display = display, width = 170, height = 30, text = "player vs computer", value = "pc"),
                Button(display = display, width = 185, height = 30, text = "computer vs computer", value = "cc")
                ],
            selectedIndex = 1
            )
    start = Button(display = display, top = 380, left = size[1]//2 - 25, width = 50, height = 30, text = "start", backgroundColor = (155,0,55))
    buttons_algorithm.drawButtons()
    difficulty.drawButtons()
    buttons_player.drawButtons()
    play_mode.drawButtons()
    start.drawButton()

    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            elif ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                if not buttons_algorithm.selectAfterCoord(pos):
                    if not difficulty.selectAfterCoord(pos):
                        if not buttons_player.selectAfterCoord(pos):
                            if not play_mode.selectAfterCoord(pos):
                                if start.selectAfterCoord(pos):
                                    display.fill((0,0,0))
                                    game.drawGameScreen()
                                    return buttons_algorithm.getValue(), difficulty.getValue(), buttons_player.getValue(), play_mode.getValue()
            pygame.display.update()

def main():
    # Initializing variables
    pygame.init()
    pygame.display.set_caption("Sociu Daniel - Quoridor")
    linesNumber = 9
    columnsNumber = 9
    cellWidth = 100
    cellHeight = 100
    cellPadding = 8
    size = (columnsNumber * (cellWidth + 1), linesNumber * (cellHeight + 1))
    screenColor = (0, 0, 0)
    screen = pygame.display.set_mode(size = size)
    background_image = pygame.image.load("background.png")
    background_image = pygame.transform.scale(background_image, size)
    screen.blit(background_image, [0, 0])

    Cell.display = screen
    Interface.initialize(linesNumber, columnsNumber, cellWidth, cellHeight, cellPadding, screenColor, screen)
    imageBlue = pygame.transform.scale(pygame.image.load("blue.png"), (Interface.imageDimension, Interface.imageDimension)) 
    imageRed = pygame.transform.scale(pygame.image.load("red.png"), (Interface.imageDimension, Interface.imageDimension)) 


    red = Piece(imageRed, 10, "red", Interface.imageDimension, Interface.linesNumber - 1, Interface.columnsNumber // 2, 0)
    blue = Piece(imageBlue, 10, "blue", Interface.imageDimension, 0, Interface.columnsNumber // 2, Interface.linesNumber - 1)
    game = Game(interface = Interface(), red = red, blue = blue)
    # Draw the menu and get the data
    algorith_type, difficulty, Game.GMIN, game_mode = draw_menu(screen, game, size)
    Game.GMAX = 'blue' if Game.GMIN == 'red' else 'red'
    MAX_DEPTH = difficulty

    current_state = State(game, "red", MAX_DEPTH)
    nextCells = []
    maxTimeComputer = 0
    minTimeComputer = 100000
    totalTimeComputer = 0
    movesComputer = 0
    start_time = None
    fontObject = pygame.font.SysFont('arial', 20)

    if game_mode == "pp":
        while True:
            if (current_state.currentPlayer == game.GMIN):
                if not start_time:
                    # Printing details / updateing screen
                    current_state.game.drawGameScreen()
                    display_turn = fontObject.render("To move: {}".format(current_state.currentPlayer), True, (255,0,0))
                    screen.blit(display_turn, display_turn.get_rect())
                    pygame.display.update()
                    start_time = int(round(time.time() * 1000))
                for ev in pygame.event.get():
                    if ev.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    if ev.type == pygame.MOUSEBUTTONDOWN and nextCells:
                        # if the user has clicked on his piece
                        pos = pygame.mouse.get_pos()
                        for cell in nextCells:
                            if cell.rectangle.collidepoint(pos):
                                user = current_state.game.pieces[current_state.currentPlayer]
                                user.line = cell.line
                                user.column = cell.column
                                nextCells = []
                                current_state.game.drawGameScreen()
                                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
                                workTime = (int(round(time.time() * 1000)) - start_time)
                                print ("User moved in: " + str(workTime))
                                start_time = None
                                if user.line == user.goal:
                                    print(user.name + " won")
                                    display_turn = fontObject.render("{} WON".format(current_state.currentPlayer), True, (255,0,0))
                                    screen.blit(display_turn, display_turn.get_rect())
                                    pygame.display.update()
                                    print("Computer minim: {}".format(minTimeComputer))
                                    print("Computer max: {}".format(maxTimeComputer))
                                    print("Computer total: {}".format(totalTimeComputer))
                                    print("Computer median: {}".format(totalTimeComputer/movesComputer))
                                    return
                                break
                        nextCells = []
                        current_state.game.drawGameScreen()
                        screen.blit(display_turn, display_turn.get_rect())
                        pygame.display.update()
                    elif ev.type == pygame.MOUSEBUTTONDOWN:
                        # Waiting for a user move, either be his piece or a wall
                        finished = False
                        pos = pygame.mouse.get_pos()
                        wallFound = []
                        currentPlayer = current_state.game.pieces[current_state.currentPlayer]
                        oppositePlayer = current_state.game.pieces[Game.GMAX]
                        if currentPlayer.moves > 0:
                            for i, line in enumerate(current_state.game.interface.cellMatrix):
                                for j, cell in enumerate(line):
                                    for k, wall in enumerate(cell.wall):
                                        if wall  and wall.collidepoint(pos):
                                            wallFound.append((cell, k, wall, i, j))
                            affectedCells = []
                            if len (wallFound) == 2:
                                affectedCells = current_state.game.getWallContinuation(wallFound)
                            if len(affectedCells) == 4 and not current_state.game.checkWallBlock(affectedCells):
                                finished = True
                                currentPlayer.moves -= 1
                                for (cell, index, wall, i, j) in affectedCells:
                                    pygame.draw.rect(current_state.game.interface.screen, cell.lineColor, wall)
                                    cell.code |= 2 ** index
                            if finished:
                                current_state.game.drawGameScreen()
                                screen.blit(display_turn, display_turn.get_rect())
                                pygame.display.update()
                                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
                                workTime = (int(round(time.time() * 1000)) - start_time)
                                print ("User moved in: " + str(workTime))
                                start_time = None
                                continue

                        user = current_state.game.pieces[current_state.currentPlayer]
                        userCell = current_state.game.interface.cellMatrix[user.line][user.column]
                        if not finished and userCell.rectangle.collidepoint(pos):
                            nextCells = current_state.game.getNextCells(current_state.currentPlayer)
                            if nextCells:
                                current_state.game.drawGameScreen(highlightCells = nextCells)
                                screen.blit(display_turn, display_turn.get_rect())
                                pygame.display.update()
            else:
                if not start_time:
                    current_state.game.drawGameScreen()
                    display_turn = fontObject.render("To move: {}".format(current_state.currentPlayer), True, (255,0,0))
                    screen.blit(display_turn, display_turn.get_rect())
                    pygame.display.update()
                    start_time = int(round(time.time() * 1000))
                for ev in pygame.event.get():
                    if ev.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    if ev.type == pygame.MOUSEBUTTONDOWN and nextCells:
                        pos = pygame.mouse.get_pos()
                        for cell in nextCells:
                            if cell.rectangle.collidepoint(pos):
                                user = current_state.game.pieces[current_state.currentPlayer]
                                user.line = cell.line
                                user.column = cell.column
                                nextCells = []
                                current_state.game.drawGameScreen()
                                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
                                workTime = (int(round(time.time() * 1000)) - start_time)
                                print ("User moved in: " + str(workTime))
                                start_time = None
                                if user.line == user.goal:
                                    print(user.name + " won")
                                    display_turn = fontObject.render("{} WON".format(current_state.currentPlayer), True, (255,0,0))
                                    screen.blit(display_turn, display_turn.get_rect())
                                    pygame.display.update()
                                    print("Computer minim: {}".format(minTimeComputer))
                                    print("Computer max: {}".format(maxTimeComputer))
                                    print("Computer total: {}".format(totalTimeComputer))
                                    print("Computer median: {}".format(totalTimeComputer/movesComputer))
                                    return
                                break
                        nextCells = []
                        current_state.game.drawGameScreen()
                        screen.blit(display_turn, display_turn.get_rect())
                        pygame.display.update()
                    elif ev.type == pygame.MOUSEBUTTONDOWN:
                        finished = False
                        pos = pygame.mouse.get_pos()
                        wallFound = []
                        currentPlayer = current_state.game.pieces[current_state.currentPlayer]
                        oppositePlayer = current_state.game.pieces[Game.GMAX]
                        if currentPlayer.moves > 0:
                            for i, line in enumerate(current_state.game.interface.cellMatrix):
                                for j, cell in enumerate(line):
                                    for k, wall in enumerate(cell.wall):
                                        if wall  and wall.collidepoint(pos):
                                            wallFound.append((cell, k, wall, i, j))
                            affectedCells = []
                            if len (wallFound) == 2:
                                affectedCells = current_state.game.getWallContinuation(wallFound)
                            if len(affectedCells) == 4 and not current_state.game.checkWallBlock(affectedCells):
                                finished = True
                                currentPlayer.moves -= 1
                                for (cell, index, wall, i, j) in affectedCells:
                                    pygame.draw.rect(current_state.game.interface.screen, cell.lineColor, wall)
                                    cell.code |= 2 ** index
                            if finished:
                                current_state.game.drawGameScreen()
                                screen.blit(display_turn, display_turn.get_rect())
                                pygame.display.update()
                                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
                                workTime = (int(round(time.time() * 1000)) - start_time)
                                print ("User moved in: " + str(workTime))
                                start_time = None
                                continue

                        user = current_state.game.pieces[current_state.currentPlayer]
                        userCell = current_state.game.interface.cellMatrix[user.line][user.column]
                        if not finished and userCell.rectangle.collidepoint(pos):
                            nextCells = current_state.game.getNextCells(current_state.currentPlayer)
                            if nextCells:
                                current_state.game.drawGameScreen(highlightCells = nextCells)
                                screen.blit(display_turn, display_turn.get_rect())
                                pygame.display.update()

    elif game_mode == "pc":
        while True:
            if (current_state.currentPlayer == game.GMIN):
                if not start_time:
                    current_state.game.drawGameScreen()
                    display_turn = fontObject.render("To move: {}".format(current_state.currentPlayer), True, (255,0,0))
                    screen.blit(display_turn, display_turn.get_rect())
                    pygame.display.update()
                    start_time = int(round(time.time() * 1000))
                for ev in pygame.event.get():
                    if ev.type == pygame.QUIT:
                        pygame.quit()
                        sys.exit()
                    if ev.type == pygame.MOUSEBUTTONDOWN and nextCells:
                        pos = pygame.mouse.get_pos()
                        for cell in nextCells:
                            if cell.rectangle.collidepoint(pos):
                                user = current_state.game.pieces[current_state.currentPlayer]
                                user.line = cell.line
                                user.column = cell.column
                                nextCells = []
                                current_state.game.drawGameScreen()
                                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
                                workTime = (int(round(time.time() * 1000)) - start_time)
                                print ("User moved in: " + str(workTime))
                                start_time = None
                                if user.line == user.goal:
                                    print(user.name + " won")
                                    display_turn = fontObject.render("{} WON".format(current_state.currentPlayer), True, (255,0,0))
                                    screen.blit(display_turn, display_turn.get_rect())
                                    pygame.display.update()
                                    print("Computer minim: {}".format(minTimeComputer))
                                    print("Computer max: {}".format(maxTimeComputer))
                                    print("Computer total: {}".format(totalTimeComputer))
                                    print("Computer median: {}".format(totalTimeComputer/movesComputer))
                                    return
                                break
                        nextCells = []
                        current_state.game.drawGameScreen()
                        screen.blit(display_turn, display_turn.get_rect())
                        pygame.display.update()
                    elif ev.type == pygame.MOUSEBUTTONDOWN:
                        finished = False
                        pos = pygame.mouse.get_pos()
                        wallFound = []
                        currentPlayer = current_state.game.pieces[current_state.currentPlayer]
                        if currentPlayer.moves > 0:
                            for i, line in enumerate(current_state.game.interface.cellMatrix):
                                for j, cell in enumerate(line):
                                    for k, wall in enumerate(cell.wall):
                                        if wall  and wall.collidepoint(pos):
                                            wallFound.append((cell, k, wall, i, j))
                            affectedCells = []
                            if len (wallFound) == 2:
                                affectedCells = current_state.game.getWallContinuation(wallFound)
                            if len(affectedCells) == 4 and not current_state.game.checkWallBlock(affectedCells):
                                finished = True
                                currentPlayer.moves -= 1
                                for (cell, index, wall, i, j) in affectedCells:
                                    pygame.draw.rect(current_state.game.interface.screen, cell.lineColor, wall)
                                    cell.code |= 2 ** index
                            if finished:
                                current_state.game.drawGameScreen()
                                screen.blit(display_turn, display_turn.get_rect())
                                pygame.display.update()
                                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
                                workTime = (int(round(time.time() * 1000)) - start_time)
                                print ("User moved in: " + str(workTime))
                                start_time = None
                                continue

                        user = current_state.game.pieces[current_state.currentPlayer]
                        userCell = current_state.game.interface.cellMatrix[user.line][user.column]
                        if not finished and userCell.rectangle.collidepoint(pos):
                            nextCells = current_state.game.getNextCells(current_state.currentPlayer)
                            if nextCells:
                                current_state.game.drawGameScreen(highlightCells = nextCells)
                                screen.blit(display_turn, display_turn.get_rect())
                                pygame.display.update()
            else:
                display_turn = fontObject.render("To move: {}".format(current_state.currentPlayer), True, (255,0,0))
                current_state.game.drawGameScreen()
                screen.blit(display_turn, display_turn.get_rect())
                pygame.display.update()
                start_time = int(round(time.time() * 1000))
                if algorith_type == "minimax":
                    new_state = min_max(current_state)
                else:
                    new_state = alpha_beta(-100, 100, current_state)
                current_state.game.blue.line = new_state.chosenState.game.blue.line
                current_state.game.blue.column = new_state.chosenState.game.blue.column
                current_state.game.blue.moves = new_state.chosenState.game.blue.moves
                current_state.game.red.line = new_state.chosenState.game.red.line
                current_state.game.red.column = new_state.chosenState.game.red.column
                current_state.game.red.moves = new_state.chosenState.game.red.moves
                for i, line in enumerate (new_state.chosenState.game.interface.cellMatrix):
                    current_line = current_state.game.interface.cellMatrix[i]
                    for j, cell in enumerate(line):
                        if (current_line[j].code != cell.code):
                            aux = byte_to_power(cell.code - current_line[j].code)
                            pygame.draw.rect(
                                current_state.game.interface.screen,
                                current_line[j].lineColor,
                                current_line[j].wall[aux]
                                )
                            current_line[j].code = cell.code
                current_state.game.drawGameScreen()
                workTime = (int(round(time.time() * 1000)) - start_time)
                print ("Computer moved in: " + str(workTime))
                maxTimeComputer = max(workTime, maxTimeComputer)
                minTimeComputer = min(workTime, maxTimeComputer)
                movesComputer += 1
                totalTimeComputer += workTime
                start_time = None
                if (current_state.game.final()):
                    print(Game.GMAX + " won")
                    display_turn = fontObject.render("{} WON".format(current_state.currentPlayer), True, (255,0,0))
                    screen.blit(display_turn, display_turn.get_rect())
                    pygame.display.update()
                    print("Computer minim: {}".format(minTimeComputer))
                    print("Computer max: {}".format(maxTimeComputer))
                    print("Computer total: {}".format(totalTimeComputer))
                    print("Computer median: {}".format(totalTimeComputer/movesComputer))
                    return
                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
    else:
        while True:
            if (current_state.currentPlayer == game.GMIN):
                display_turn = fontObject.render("To move: {}".format(current_state.currentPlayer), True, (255,0,0))
                current_state.game.drawGameScreen()
                screen.blit(display_turn, display_turn.get_rect())
                pygame.display.update()
                start_time = int(round(time.time() * 1000))
                if algorith_type == "minimax":
                    new_state = min_max(current_state)
                else:
                    new_state = alpha_beta(-100, 100, current_state)
                current_state.game.blue.line = new_state.chosenState.game.blue.line
                current_state.game.blue.column = new_state.chosenState.game.blue.column
                current_state.game.blue.moves = new_state.chosenState.game.blue.moves
                current_state.game.red.line = new_state.chosenState.game.red.line
                current_state.game.red.column = new_state.chosenState.game.red.column
                current_state.game.red.moves = new_state.chosenState.game.red.moves
                for i, line in enumerate (new_state.chosenState.game.interface.cellMatrix):
                    current_line = current_state.game.interface.cellMatrix[i]
                    for j, cell in enumerate(line):
                        if (current_line[j].code != cell.code):
                            aux = byte_to_power(cell.code - current_line[j].code)
                            pygame.draw.rect(
                                current_state.game.interface.screen,
                                current_line[j].lineColor,
                                current_line[j].wall[aux]
                                )
                            current_line[j].code = cell.code
                workTime = (int(round(time.time() * 1000)) - start_time)
                if workTime < 100:
                    time.sleep(0.3)
                current_state.game.drawGameScreen()
                print ("Computer moved in: " + str(workTime))
                maxTimeComputer = max(workTime, maxTimeComputer)
                minTimeComputer = min(workTime, maxTimeComputer)
                movesComputer += 1
                totalTimeComputer += workTime
                start_time = None
                if (current_state.game.final()):
                    print(Game.GMAX + " won")
                    display_turn = fontObject.render("{} WON".format(current_state.currentPlayer), True, (255,0,0))
                    screen.blit(display_turn, display_turn.get_rect())
                    pygame.display.update()
                    print("Computer minim: {}".format(minTimeComputer))
                    print("Computer max: {}".format(maxTimeComputer))
                    print("Computer total: {}".format(totalTimeComputer))
                    print("Computer median: {}".format(totalTimeComputer/movesComputer))
                    return
                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)
            else:
                display_turn = fontObject.render("To move: {}".format(current_state.currentPlayer), True, (255,0,0))
                current_state.game.drawGameScreen()
                screen.blit(display_turn, display_turn.get_rect())
                pygame.display.update()
                start_time = int(round(time.time() * 1000))
                if algorith_type == "minimax":
                    new_state = min_max(current_state, 2)
                else:
                    new_state = alpha_beta(-100, 100, current_state, 2)
                current_state.game.blue.line = new_state.chosenState.game.blue.line
                current_state.game.blue.column = new_state.chosenState.game.blue.column
                current_state.game.blue.moves = new_state.chosenState.game.blue.moves
                current_state.game.red.line = new_state.chosenState.game.red.line
                current_state.game.red.column = new_state.chosenState.game.red.column
                current_state.game.red.moves = new_state.chosenState.game.red.moves
                for i, line in enumerate (new_state.chosenState.game.interface.cellMatrix):
                    current_line = current_state.game.interface.cellMatrix[i]
                    for j, cell in enumerate(line):
                        if (current_line[j].code != cell.code):
                            aux = byte_to_power(cell.code - current_line[j].code)
                            pygame.draw.rect(
                                current_state.game.interface.screen,
                                current_line[j].lineColor,
                                current_line[j].wall[aux]
                                )
                            current_line[j].code = cell.code
                workTime = (int(round(time.time() * 1000)) - start_time)
                if workTime < 100:
                    time.sleep(0.3)
                current_state.game.drawGameScreen()
                print ("Computer moved in: " + str(workTime))
                maxTimeComputer = max(workTime, maxTimeComputer)
                minTimeComputer = min(workTime, maxTimeComputer)
                movesComputer += 1
                totalTimeComputer += workTime
                start_time = None
                if (current_state.game.final()):
                    print(Game.GMAX + " won")
                    display_turn = fontObject.render("{} WON".format(current_state.currentPlayer), True, (255,0,0))
                    screen.blit(display_turn, display_turn.get_rect())
                    pygame.display.update()
                    print("Computer minim: {}".format(minTimeComputer))
                    print("Computer max: {}".format(maxTimeComputer))
                    print("Computer total: {}".format(totalTimeComputer))
                    print("Computer median: {}".format(totalTimeComputer/movesComputer))
                    return
                current_state.currentPlayer = current_state.game.getOppositePlayer(current_state.currentPlayer)

if __name__ == "__main__":
    main()
    while True:
        for ev in pygame.event.get():
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()

