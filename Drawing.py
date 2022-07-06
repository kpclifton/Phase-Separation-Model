import sys
import time
import pandas as pd
import random
import numpy as np
from pylab import *
from scipy.ndimage import measurements
from MoveProtein import MoveProtein
from PyQt5.QtGui import QPainter, QColor, QFont, QPen, QBrush, QPolygonF
from PyQt5.QtCore import Qt, QRect, QLine, QPoint, QRectF, QCoreApplication
from PyQt5.QtWidgets import QWidget, QApplication, QPushButton, QInputDialog, QLineEdit,QComboBox

# Size of windows in pixels is = size + AxisSpace (offset)
windowSizeX = 1000
windowSizeY = 500
yAxisSpace = 100
xAxisSpace = 40



class AggregationModel(QWidget):

  def __init__(self):
    super().__init__()
    #GFX and constant inits
    self.setGeometry(300,50, windowSizeX+xAxisSpace, windowSizeY + yAxisSpace )
    self.setWindowTitle('Aggregation Model')

    #Adjustable variables
    self.destabilizationConstant = 2 #How much to change temperature in destabilized area (multiplicative)
    self.record = False #Whether or not we should record clustering over time


    #Constants
    self.gridSizeX = 100 #Number of squares in width. Grid is 2x as wide as is tall.
    self.gridSizeY = self.gridSizeX//2 #Rectangular cell

    #Flags
    self.destabilizationActive = False
    self.simulating = False #Are we currently simulating

    #Set up function calls
    self.gridSetup()
    self.dataInit()
    self.uiInit()
    self.show()

  def dataInit(self):
    #Main Data Arrays
    initialProtein = 700
    temp = ([0]*((self.gridSizeX*self.gridSizeY)-initialProtein))+([1]*initialProtein)
    self.occupationData = np.random.choice(temp, size = (self.gridSizeY,self.gridSizeX), replace = False) #Randomly generate the locations of 700 proteins
    self.temperatureData = np.ones((self.gridSizeY,self.gridSizeX)) #array of ones
    #self.temperatureData = np.ones((self.gridSizeY,self.gridSizeX))*2 #array for different initial temperatures
    self.clusterAvg = []
    self.clusterAvgleft = []
    self.clusterAvgright = []
    del temp




  def destabilization(self):
      if self.destabilizationActive:
          self.destabilizationActive = False
          self.destabilizationBtn.setText('Destabilize')
          self.temperatureData = np.ones((self.gridSizeY,self.gridSizeX)) #array of ones
          #self.temperatureData = np.ones((self.gridSizeY,self.gridSizeX))*2 #array of twos
      else:
          self.destabilizationActive = True
          self.destabilizationBtn.setText('Stabilize')
          #This is where you can input alternative destabilization regimes
          #In current form destabilizes half of the cell
          for i in range(self.gridSizeY): #Yes I know this is inefficient, but it does not matter as we do it one time. Would be best to use either mapping or at least do the numpy iter call.
            self.temperatureData[i,(self.gridSizeX//2):] = 1*self.destabilizationConstant



  def calculateClusters(self):
      #calculate mean clusters of total
      lw, num = measurements.label(self.occupationData,[[1,1,1],[1,1,1], [1,1,1]])
      area = measurements.sum(self.occupationData, lw, index=arange(1,lw.max() + 1))
      self.clusterAvg.append(round(np.average(area[1:]),3))

      #calculate mean clusters of left side
      lw2, num2 = measurements.label(self.occupationData[:,:(self.gridSizeX//2)], [[1,1,1],[1,1,1], [1,1,1]])
      area2 = measurements.sum(self.occupationData[:,:(self.gridSizeX//2)], lw2, index=arange(1,lw2.max() + 1))
      self.clusterAvgleft.append(round(np.average(area2[1:]),3))

      #calculate mean clusters of right side
      lw3, num3 = measurements.label(self.occupationData[:,(self.gridSizeX//2):], [[1,1,1],[1,1,1], [1,1,1]])
      area3 = measurements.sum(self.occupationData[:,(self.gridSizeX//2):], lw3, index=arange(1,lw3.max() + 1))
      self.clusterAvgright.append(round(np.average(area3[1:]),3))
#After this point mostly backend stuff
  def paintEvent(self, event):
    qp = QPainter()
    blackPen = QPen(QBrush(Qt.black),0)
    qp.begin(self)
    qp.setPen(blackPen)
    #Draws Steps Taken
    qp.drawText(750,30,"Steps taken: "+str(self.stepsTaken))
    #Fills in squares with white if they are filled, otherwise with black.
    for row in range(self.gridSizeY):
        for column in range(self.gridSizeX):
            if self.occupationData[row][column]:
                qp.fillRect(self.board[column][row],Qt.white)
            else:
                qp.fillRect(self.board[column][row],Qt.black)

    #Hack for getting sequential updating
    if self.batchNumber >=0 and self.simulating ==True:
        self.tick()
        if self.batchNumber <= 0:
            self.simulating = False
            self.goBtn.setEnabled(True)
        self.batchNumberButton.setText(str(self.batchNumber))


    qp.end()



  def saveCSV(self):
      pd.DataFrame(self.clusterAvg).to_csv("clusteringAvgs")
      pd.DataFrame(self.clusterAvgleft).to_csv("clusteringAvgsleft")
      pd.DataFrame(self.clusterAvgright).to_csv("clusteringAvgsright")

  def gridSetup(self):
    #Constants
    self.board = []
    #sets up the graphics list, as well as the size of squares
    self.squareSize = windowSizeX/self.gridSizeX
    self.board = [[0 for row in range(self.gridSizeY)] for column in range(self.gridSizeX)]
    for row in range(0,self.gridSizeY):
      for column in range(0,self.gridSizeX):
        self.board[column][row] = QRect(column * self.squareSize + xAxisSpace/2,row*self.squareSize + yAxisSpace/2,self.squareSize,self.squareSize)


#UI stuff

  def uiInit(self):
    self.batchNumber = 10 #Number of batches performed when simulation is run
    self.batchSize = 1 #Number of steps to do in a given batch (before updating screen)
    self.stepsTaken = 0

    self.goBtn = QPushButton("Go", self)
    self.goBtn.clicked.connect(self.timekeeper)
    self.goBtn.resize(30,30)
    self.goBtn.move(30,10)

    self.destabilizationBtn = QPushButton("Destabilize", self)
    self.destabilizationBtn.clicked.connect(self.destabilization)
    self.destabilizationBtn.resize(120,30)
    self.destabilizationBtn.move(100,10)

    self.batchNumberButton = QPushButton(str(self.batchNumber), self)
    self.batchNumberButton.clicked.connect(self.getBatchNumberInteger)
    self.batchNumberButton.resize(120,30)
    self.batchNumberButton.move(250,10)

    self.batchSizeButton = QPushButton(str(self.batchSize), self)
    self.batchSizeButton.clicked.connect(self.getBatchSizeInteger)
    self.batchSizeButton.resize(80,30)
    self.batchSizeButton.move(400,10)
    if self.record:
        self.saveBtn = QPushButton("Save", self)
        self.saveBtn.clicked.connect(self.saveCSV)
        self.saveBtn.resize(100,30)
        self.saveBtn.move(600,10)

  def getBatchNumberInteger(self):
        i, okPressed = QInputDialog.getInt(self, "Enter an intger","Batches to do", 1, 0, 10000000, 1)
        if okPressed:
            self.batchNumber = i
            self.batchNumberButton.setText(str(i))

  def getBatchSizeInteger(self):
        i, okPressed = QInputDialog.getInt(self, "Enter an intger","Batch Size", 1, 0, 10000, 1)
        if okPressed:
            self.batchSize = i
            self.batchSizeButton.setText(str(i))
  def timekeeper(self):
       if self.batchNumber<=0:
           self.goBtn.setEnabled(True)
       else:
           self.simulating = True
           self.goBtn.setEnabled(False)
           self.update()

  def tick(self):
      if self.batchNumber <=0:
          return
      for i in range(self.batchSize):
          MoveProtein(self.occupationData,self.temperatureData,self.gridSizeX,self.gridSizeY)
          if self.record:
            self.calculateClusters()
          self.stepsTaken+=1
      self.batchNumber -=1
      self.update()


if __name__ == '__main__':
  app = QApplication(sys.argv)
  ex = AggregationModel()
  sys.exit(app.exec_())
