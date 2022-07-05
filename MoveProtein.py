import random
import numpy as np
def MoveProtein(occupationData,temperatureData,gridSizeX,gridSizeY):
  # step 1: pick random occupied position
  position = random.choice(np.argwhere(occupationData == 1))

  # save indices of randomly picked occupied position (i,j)
  i = position[0]
  j = position[1]

  #step 2: pick random potential direction for diffusion
  #(while obeying the reflective boundary)

  #initialize list with all possible directions
  direction = ["A","B","C","D","E","F","G","H"]

  #remove directions that are not possible based on if the position of interest is
  #a middle, edge, or corner spot

  # bounds of grid
  I = gridSizeY
  J = gridSizeX

  if ((i - 1) < 0) or ((j-1) < 0):
   direction.remove("A")

  if (i - 1) < 0:
   direction.remove("B")

  if ((i-1) <0) or ((j+1) == J):
   direction.remove("C")

  if ((j+1) == J):
   direction.remove("D")

  if ((i+1) == I) or ((j+1) == J):
   direction.remove("E")

  if ((i+1) == I):
   direction.remove("F")

  if ((i+1) == I) or ((j-1) < 0):
   direction.remove("G")

  if ((j-1) < 0):
   direction.remove("H")

  #choose a random direction from the directions that are possible
  r_2 = random.choice((direction))

  #step 3: check occupancy of square in that direction and save the position of that
  #neighboring square

  if r_2 == "A":
    occupancy = occupationData[i-1][j-1]
    direction_to_move = [i-1,j-1]

  if r_2 == "B":
    occupancy = occupationData[i-1][j]
    direction_to_move = [i-1,j]

  if r_2 == "C":
    occupancy = occupationData[i-1][j+1]
    direction_to_move = [i-1,j+1]

  if r_2 == "D":
    occupancy = occupationData[i][j+1]
    direction_to_move = [i,j+1]

  if r_2 == "E":
    occupancy = occupationData[i+1][j+1]
    direction_to_move = [i+1,j+1]

  if r_2 == "F":
    occupancy = occupationData[i+1][j]
    direction_to_move = [i+1,j]

  if r_2 == "G":
    occupancy = occupationData[i+1][j-1]
    direction_to_move = [i+1,j-1]

  if r_2 == "H":
    occupancy = occupationData[i][j-1]
    direction_to_move = [i,j-1]

  #step 4: determine if diffuse in that direction
  #if unoccupied at neighboring square in the direction to potentially move, count how many
  #occupied neighbors there are for position (i,j). This count gives number of bonds that
  #would be lost if the monomer at position (i,j) diffused
  if occupancy == 0:
    n_lost = 0
    for p in range(0,len(direction)):
      if direction[p] == "A":
        n_lost = n_lost + occupationData[i-1][j-1]
      if direction[p] == "B":
        n_lost = n_lost + occupationData[i-1][j]
      if direction[p] == "C":
        n_lost = n_lost + occupationData[i-1][j+1]
      if direction[p] == "D":
        n_lost = n_lost + occupationData[i][j+1]
      if direction[p] == "E":
        n_lost = n_lost + occupationData[i+1][j+1]
      if direction[p] == "F":
        n_lost = n_lost + occupationData[i+1][j]
      if direction[p] == "G":
        n_lost = n_lost + occupationData[i+1][j-1]
      if direction[p] == "H":
        n_lost = n_lost + occupationData[i][j-1]

    #calculate reaction rate for diffusion and unbinding
    k_0 = 1
    E = 1
    temp = temperatureData[i][j]

    k = k_0*np.exp((-1)*E*n_lost/temp)

    #generate random number in [0,1)
    r_3 = np.random.rand()

    # if reaction rate is greater than the random number then diffuse in direction chosen
    if k > r_3:
      occupationData[i][j] = 0
      occupationData[direction_to_move[0]][direction_to_move[1]] = 1
