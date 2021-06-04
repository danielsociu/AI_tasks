import time
import copy
import pygame
import sys

class Cell:
    wallWidth = 15
    backgroundCell = (255,255,255)
    lineColor = (33, 222, 71)

    def __init__(self, display, interface,  left, top, width, height, line, column, code = 0):
        self.rectangle = pygame.Rect(left, top, width, height)
        self.display = display
        self.wall = [None, None, None, None]
        self.code = 0
        if line > 0:
            self.wall[0] = pygame.Rect(left, top - 1 - self.__class__.wallWidth // 2,\
                    width, self.__class__.wallWidth)
        else:
            self.code += 2**0
        if column < interface.columnsNumber - 1:
            self.wall[1] = pygame.Rect(left + width - self.__class__.wallWidth // 2, top,\
                    self.__class__.wallWidth, height)
        else:
            code += 2**1
        if line < interface.linesNumber - 1:
            self.wall[2] = pygame.Rect(left, top + height - self.__class__.wallWidth // 2, \
                    width, self.__class__.wallWidth)
        else:
            code += 2**2
        if column > 0:
            self.wall[3] = pygame.Rect(left - self.__class__.wallWidth // 2, top, \
                    self.__class__.wallWidth, height)
        else:
            code += 2**3

    def drawCell(self):
        pygame.draw.rect(self.display, self.__class__.backgroundCell, self.rectangle)
        byte = 1
        for i in range(4):
            if self.code & byte:
                if self.wall[i]:
                    pygame.draw.rect(self.display, self.__class__.lineColor, self.wall[i])
            byte *= 2

class Piece:
    def __init__(self, path, dimension, line, column):
        self.path = path 
        self.dimension = dimension
        self.line = line
        self.column = column
        image = pygame.image.load(path)
        self.image = pygame.transform.scale(image, (dimension, dimension))


class Interface:
    screenColor = (0, 0, 0)

    def __init__(self, linesNumber, columnsNumber, cellWidth, cellHeight, cellPadding, screen = None):
        self.linesNumber = linesNumber
        self.columnsNumber = columnsNumber
        self.cellWidth = cellWidth
        self.cellHeight = cellHeight
        self.cellPadding = cellPadding
        self.imageDimension = min(cellWidth, cellHeight) - 2 * cellPadding
        self.screen = screen
        self.cellMatrix = [[
            Cell(
                screen,
                self, 
                col * (self.cellWidth + 1),
                line * (self.cellHeight + 1),
                self.cellWidth,
                self.cellHeight,
                line,
                col
                ) for col in range(columnsNumber)]
            for line in range(linesNumber)
            ]

    def drawImage(self, image, cell):
        self.screen.blit(
                image,
                (
                    cell.rectangle.left + self.cellPadding,
                    cell.rectangle.top + self.cellPadding
                ))

    def drawGameScreen(self, pieces):
        self.screen.fill(self.screenColor)
        for i, line in enumerate(self.cellMatrix):
            for j, cell in enumerate(line):
                cell.drawCell()
                for piece in pieces:
                    if i == piece.line and j == piece.column:
                        self.drawImage(piece.image, cell)
        pygame.display.update()


class Game:
    GMIN = None
    GMAX = None

    def __init__(self, interface):
        self.red = Piece("red.png", interface.imageDimension, interface.linesNumber - 1, interface.columnsNumber // 2)
        self.blue = Piece("blue.png", interface.imageDimension, 0, interface.columnsNumber // 2)
        self.pieces = [self.red, self.blue]
        self.interface = interface

    def drawGameScreen(self):
        self.interface.drawGameScreen(self.pieces)


    def getOppositeWalls(self, x, y, index):
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
        nextIndex1 = self.getOppositeWalls(x, y, index)
        nextIndex2 = self.getOppositeWalls(x, y, (index + 2) % 4)
        if cell1.code & 2**nextIndex1 and cell2.code & 2**nextIndex2:
            return False
        return True


    def checkValidationWall(self, i, j, moveX, moveY, index):
        nextCell = self.interface.cellMatrix[i + moveX][j + moveY]
        if (nextCell.code & 2**index):
            return None
        return (nextCell, index, nextCell.wall[index], i + moveX, j + moveY)

    def getWallContinuation(self, currentWall):
        affectedCells = []
        bad = 0
        if (currentWall[0][1] & 1) != (currentWall[1][1] & 1) or (currentWall[0][0].code & 2**currentWall[0][1]):
            return currentWall
        for other, (cell, index, wall, i, j) in enumerate(currentWall):
            nextTuple = None
            otherCell = currentWall[other ^ 1][0]
            if index & 1:
                if i < self.interface.linesNumber - 1 and self.checkOppositeWall(cell, otherCell, 1, 0, index):
                    move = 1
                    otherTuple = currentWall[other ^ 1]
                    nextTuple = self.checkValidationWall(i, j, move, 0, index)
                if not nextTuple and i > 0 and self.checkOppositeWall(cell, otherCell, -1, 0, index):
                    move = -1
                    nextTuple = self.checkValidationWall(i, j, move, 0, index)
            else:
                if j < self.interface.columnsNumber - 1 and self.checkOppositeWall(cell, otherCell, 0, 1, index):
                    move = 1
                    nextTuple = self.checkValidationWall(i, j, 0, move, index)
                if not nextTuple and j > 0 and self.checkOppositeWall(cell, otherCell, 0, -1, index):
                    move = -1
                    nextTuple = self.checkValidationWall(i, j, 0, move, index)

            if nextTuple != None:
                affectedCells.append(nextTuple)
        affectedCells += currentWall
        return affectedCells

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
    buttons_player = ButtonGroup(
            top = 200,
            left = size[1]//2 - 60,
            buttonList = [
                Button(display = display, width = 50, height = 30, text = "red", value = "red"),
                Button(display = display, width = 50, height = 30, text = "blue", value = "blue")
                ]
            )
    play_mode = ButtonGroup(
            top = 260,
            left = size[1]//2 - (190 * 3) / 2,
            buttonList = [
                Button(display = display, width = 160, height = 30, text = "player vs player", value = "pp"),
                Button(display = display, width = 170, height = 30, text = "player vs computer", value = "pc"),
                Button(display = display, width = 185, height = 30, text = "computer vs computer", value = "cc")
                ],
            selectedIndex = 1
            )
    start = Button(display = display, top = 320, left = size[1]//2 - 25, width = 50, height = 30, text = "start", backgroundColor = (155,0,55))
    buttons_algorithm.drawButtons()
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
                    if not buttons_player.selectAfterCoord(pos):
                        if not play_mode.selectAfterCoord(pos):
                            if start.selectAfterCoord(pos):
                                display.fill((0,0,0))
                                game.drawGameScreen()
                                return buttons_algorithm.getValue(), buttons_player.getValue(), play_mode.getValue()
            pygame.display.update()

def main():
    pygame.init()
    pygame.display.set_caption("Quoridor")
    linesNumber = 9
    columnsNumber = 9
    cellWidth = 100
    cellHeight = 100
    cellPadding = 8
    size = (columnsNumber * (cellWidth + 1), linesNumber * (cellHeight + 1))
    screen = pygame.display.set_mode(size = size)
    background_image = pygame.image.load("background.png")
    background_image = pygame.transform.scale(background_image, size)
    screen.blit(background_image, [0, 0])
    game = Game(
            Interface(linesNumber, columnsNumber, cellWidth, cellHeight, cellPadding, screen)
            )
    algorith_type, Game.GMIN, game_mode = draw_menu(screen, game, size)

    while True:
        for ev in pygame.event.get():
            game.drawGameScreen()
            if ev.type == pygame.QUIT:
                pygame.quit()
                sys.exit()
            if ev.type == pygame.MOUSEBUTTONDOWN:
                pos = pygame.mouse.get_pos()
                wallFound = []
                for i, line in enumerate(game.interface.cellMatrix):
                    for j, cell in enumerate(line):
                        for k, wall in enumerate(cell.wall):
                            if wall  and wall.collidepoint(pos):
                                wallFound.append((cell, k, wall, i, j))
                affectedCells = []
                if len (wallFound) == 2:
                    affectedCells = game.getWallContinuation(wallFound)
                if len(affectedCells) == 4:
                    for (cell, index, wall, i, j) in affectedCells:
                        pygame.draw.rect(game.interface.screen, cell.lineColor, wall)
                        cell.code |= 2 ** index

if __name__ == "__main__":
    main()

