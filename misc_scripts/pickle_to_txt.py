import pickle
import matplotlib.pyplot as plt

with open("triangle_atom_assignment_upper_test2.pickle", "rb") as f:
     data = pickle.load(f)

#turn pickle into a text file, so I can visualize the data
with open("output.txt", "w") as f:
    f.write(str(data))

