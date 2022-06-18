import os
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

os.chdir("C:\\Users\\antony\\Desktop\\suhita lab\\hodgkin-huxley\\run_data\\test")
path = os.getcwd()
file = pd.read_csv("filename.csv")
print(file)
csv_files = file["file names"]
# csv_files.sort_values()
print(csv_files)
data = []
for file in csv_files:
    data.append(pd.read_csv(file))

# print(data[80]['Membrane Voltage (mV)'])
# print(data[80]['Time (ms)'])

max_voltages = np.zeros(300, float)
min_voltages = np.zeros(300, float)
i = 0

for frame in data:
    # fig, ax = plt.subplots()
    max_voltages[i] = frame['Membrane Voltage (mV)'][-2000:].max()
    print(max_voltages[i])
    min_voltages[i] = frame['Membrane Voltage (mV)'][-2000:].min()
    # print(min_voltages[i])
    i = i + 1

input_current = np.arange(0,300, 1)

plt.plot(input_current, max_voltages, label="Maximum Voltage")
plt.plot(input_current, min_voltages, label="Minimum Voltage")
plt.xlabel('Input current (micro Amps/cm^2)')
plt.ylabel('Voltage (mV)')
plt.title('Bifurcation diagram')
plt.grid()
plt.legend()
plt.show()

# df = pd.DataFrame({"Input Current": input_current, "max Membrane Voltage (mV)": max_voltages, "min Membrane Voltage ("
#                                                                                               "mV)": min_voltages})
# df.to_csv("output_max and min voltages.csv", index=False)

# # fig, ax = plt.subplots()
# ax.plot(input_current, max_voltages)
# ax.set_xlabel('Input current')
# ax.set_ylabel('Voltage (mV)')
# ax.set_title('Bifurcation diagram')
# plt.grid()
#
# # fig, ax = plt.subplots()
# ax.plot(input_current, min_voltages)
# ax.set_xlabel('Input current')
# ax.set_ylabel('Voltage (mV)')
# ax.set_title('Bifurcation diagram')
# plt.grid()



# plt.show()
