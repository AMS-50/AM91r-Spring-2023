import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
import pandas as pd

def initial_data_display():
    #read the file we want to use
    df = pd.read_csv('det_5sigma.csv')
    print(df)

    #
    plt.scatter(df.col0 - 360, df.col1, s=1, c=np.tan(df.col2))
    plt.xlim(-8.9, -9.2)
    plt.show()

    min = df.col2[0]
    max = df.col2[0]
    last_time = df.col2[0]

    time_step_counter = 0
    for time in df.col2:
        if time > max:
            max = time
        if time < min:
            min = time
        if last_time < time:
            last_time = time
            time_step_counter += 1

    print("min = ",min)
    print("max = ",max)
    print("dif = ",max - min)
    # +1 because number of steps is changes + 1
    print("num of time steps", time_step_counter + 1)



def guess_propogate():
    #read file
    df = pd.read_csv('det_5sigma.csv')
    #display initial region
    
    
    plt.scatter(df.col0 - 360, df.col1, s=1, c=np.tan(df.col2))
    plt.xlim(-8.9, -9.2)
    plt.xlim(-8.91, -8.96)
    plt.ylim(-4.99, -4.97)
    plt.show()



    # calculate first time 
    # assumes detection sorted in time order
    first_time = df.col2[0]
    # manual alpha beta guesses
    alpha = 0.12
    beta = 0.056
    plt.scatter(df.col0 - 360 + alpha * (df.col2 - first_time), df.col1 + beta * (df.col2 - first_time), s=1, c=np.tan(df.col2))
    plt.xlim(-8.9, -9.2)
    plt.xlim(-8.91, -8.96)
    plt.ylim(-4.99, -4.97)
    plt.show()
    return 


def closest_detection(theta_x, theta_y, file_name):

    df = pd.read_csv(file_name)

    x = df.col0[0] - 360
    y = df.col1[0]
    time = df.col2[0]
    dist = np.sqrt(np.square(theta_x - x) + np.square(theta_y - y))
    print(dist)

    for index, row in df.iterrows():
        new_dist = np.sqrt(np.square(theta_x - (row['col0'] - 360)) + np.square(theta_y - row['col1']))
        #print(theta_x, row['col0'] - 360 , theta_y, row['col1'], new_dist)
        if new_dist < dist:
            dist = new_dist
            x = row['col0'] - 360
            y = row['col1']
            time = row['col2']


    return (x,y,time)

def linearity_checker():
    #manual input for now
    file_name = 'det_5sigma.csv'
    df = pd.read_csv(file_name)
    #manually input for now
    head = closest_detection(-9.0740, -4.9512, file_name)
    tail = closest_detection(-9.1063, -4.9737, file_name)
    print (head)
    print (tail)


    # create a linear function from head to tail
    alpha_dot = (tail[0] - head[0]) / (tail[2] - head[2])
    beta_dot = (tail[1] - head[1]) / (tail[2] - head[2])
    # plt.xlim(head[0], tail[0])
    # plt.ylim(tail[1], head[1])
    # plt.scatter(df.col0 - 360, df.col1, s=1, c=np.tan(df.col2))
    # plt.plot(head[0], tail[0], marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")
    # plt.plot(head[1], tail[1], marker="o", markersize=20, markeredgecolor="red", markerfacecolor="green")
    # plt.show()


    # plot for the difference

    plt.ylim(0, 0.0002)
    plt.xlim(head[2], tail[2])
    plt.scatter(df.col2,
        np.sqrt
        (np.square(df.col0 - 360 - (head[0] + (alpha_dot * (df.col2 - head[2]))))
        + np.square(df.col1 - (head[1] + (beta_dot * (df.col2 - head[2]))))),
         s=1, c=np.tan(df.col2))
    plt.show()
    


    # plot points on the linear function if thy are within 0.3 of where they should be theoritcally


    return

# calculates residual of detections to some line
def residual(alpha, alpha_dot, beta, beta_dot, gammas, x, y, t):
    x_loc = alpha + alpha * t + XE(t)
    y_loc = beta + beta * t + BE(t)
    distance = np.sqrt(
        np.square(x_loc - x) 
        + np.square(y_loc - y)
    )
    #TO DO
    return

#order of opperations
#guess gamma
#guess alpha dot, beta dot
#back propogate
#looking at binnings and go to the highest (itterte over largest to least large)
#if hi

def height_and_width_check(file):
    max_height = df.col0[0]
    min_height = df.col0[0]
    max_width = df.col1[0] - 360
    min_width = df.col1[0] - 360
    df = pd.read_csv(file)
    for index, row in df.iterrows():
        # calculate new values
        if (row['col0'] - 360 < min_width):
            min_width = row['col0'] - 360
        if (row['col0'] - 360 > max_width):
            max_width = row['col0'] - 360
        if (row['col1'] < min_height):
            min_height = row['col1']
        if (row['col1'] > max_height):
            max_height = row['col1']

    return(max_width - min_width, max_height - min_height, min_width, min_height)

def manual_algo():
    file = 'det_5sigma.csv'
    # read the file we want to use
    df = pd.read_csv(file)

    # bin sizes and data set construction
    # get extremitites
    (total_theta_x, total_theta_y, start_x, start_y) = height_and_width_check(file)
    #PLACE HOLDER NUMBERS
    x_bin_size = 1
    y_bin_size = 1
    bin = np.zeroes(np.ceil(total_theta_x / x_bin_size), np.ceil(total_theta_y / y_bin_size))

    # pick our variables (will be automated later)
    alpha_dot = 0
    beta_dot = 0
    gamma = 0

    # create bins

    # itterate to give guesses
    for index, row in df.iterrows():
        # calculate new values
        new_theta_x = row['col0'] - 360 + (alpha_dot * row['col2']) #+ gamma * function of time
        new_theta_y = row['col1'] + (beta_dot * row['col2']) #+ gamma * function of time

        # increase bin count
        bin[np.floor((start_x - new_theta_x) /x_bin_size)][np.floor((start_y - new_theta_y) /y_bin_size)] += 1


    # figure out highest bin and neighbours within distance d
    d = 1
    # for now just plan to re-itterate over all of them
    # TO DO

    # run the regression on them
    # plan
    # least square fucntion scipy.optimize.least_squares
    # run on alpha, alpha dot, gamma, beta and beta dot
    # residual is given some time t' for the detection, its distance to the line's point at t'
    # x_0 is a matrix of guesses

    #EXAMPLE FOR NON LINEAR 
    #https://scipy-lectures.org/intro/summary-exercises/optimize-fit.html

    

def main():
    #initial_data_display()
    guess_propogate()
    #linearity_checker()
    return

if __name__ == "__main__":
    main()


