# Original author: Kevin
# Modified by: Binyamin
# Modified further by : Glenn - 24/03/2021
# And even further by: Eric - Summer/2021

#!/usr/bin/env python
# coding: utf-8

import itertools, threading, time, sys, os, re
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import shutil

progress = True
log_file_folder = ""

# Handles command line flags
for flag in sys.argv:
    lowered = flag.lower()
    if lowered == "--no-progress" or lowered == "-np":
        progress = False
    elif "-log-dir" in lowered or "-ld" in lowered:
        log_file_folder = flag.split("=",1)[1]

if log_file_folder == "":
    print("\n\033[31mASSERT:\033[m Must set a log folder path using the flag -log-dir=<path to logs folder>\033[0m")
    exit(-1)

# Loading animation
done    = False
success = True
def animate():
    # Loop through the animation cycles
    for c in itertools.cycle(["|", "/", "-", "\\"]):
        # When the global variable is set at the eof
        #  then the infinite loop will break
        if done:
            break
        # Using stdout makes it easier to deal with the newlines inherent in print()
        # '\r' replaces the previous line so we don't crowd the terminal
        sys.stdout.write("\r\033[33maggregating graphs " + c + "\033[0m")
        sys.stdout.flush()
        time.sleep(0.1)
    if success:
        sys.stdout.write("\r\033[1;32mDone.                   \033[0m\n")

# Don't forget to thread it!
if progress:
    t = threading.Thread(target=animate)
    t.start()

# Setup paths, filenames, and folders
log_filename = log_file_folder + "/pandemic_state.txt"
path         = log_file_folder + "/stats/aggregate"
base_name    = path + "/"
shutil.rmtree(path, ignore_errors=True)

# Regex str to find underscore and one or more characters after the underscore (model id)
regex_model_id = "_\w+"
# Regex str to read all state contents between <>
regex_state = "<.+>"

# State log structure
sIndex    = 1
eIndex    = 2
vd1Index  = 3
vd2Index  = 4
iIndex    = 5
rIndex    = 6
neweIndex = 7
newiIndex = 8
newrIndex = 9
dIndex    = 10

COLOR_SUSCEPTIBLE = 'xkcd:blue'
COLOR_INFECTED    = 'xkcd:red'
COLOR_EXPOSED     = 'xkcd:sienna'
COLOR_DOSE1       = '#B91FDE'
COLOR_DOSE2       = '#680D5A'
COLOR_RECOVERED   = 'xkcd:green'
COLOR_DEAD        = 'xkcd:black'

STYLE_SUSCEPTIBLE = '-'
STYLE_INFECTED    = '--'
STYLE_EXPOSED     = '-.'
STYLE_DOSE1       = 'solid'
STYLE_DOSE2       = 'dashed'
STYLE_RECOVERED   = ':'
STYLE_DEAD        = 'solid'

data        = []
curr_time   = None
curr_states = {}
total_pop   = {}

def curr_states_to_df_row(sim_time, curr_states, total_pop):
    total_S = total_E = total_VD1 = total_VD2 = total_I = total_R = total_D = 0
    new_E = new_I = new_R = total_S

    # Sum the number of S,E,V,I,R,D persons in all cells
    for key in curr_states:
        cell_population = curr_states[key][0]
        total_S     += round(cell_population*(curr_states[key][sIndex]))
        total_E     += round(cell_population*(curr_states[key][eIndex]))
        total_VD1   += round(cell_population*(curr_states[key][vd1Index]))
        total_VD2   += round(cell_population*(curr_states[key][vd2Index]))
        total_I     += round(cell_population*(curr_states[key][iIndex]))
        total_R     += round(cell_population*(curr_states[key][rIndex]))
        total_D     += round(cell_population*(curr_states[key][dIndex]))

        new_E += round(cell_population*(curr_states[key][neweIndex]))
        new_I += round(cell_population*(curr_states[key][newiIndex]))
        new_R += round(cell_population*(curr_states[key][newrIndex]))

    # then divide each by the total population to get the percentage of population in each state
    percent_S   = total_S   / total_pop
    percent_E   = total_E   / total_pop
    percent_VD1 = total_VD1 / total_pop
    percent_VD2 = total_VD2 / total_pop
    percent_I   = total_I   / total_pop
    percent_R   = total_R   / total_pop
    percent_D   = total_D   / total_pop

    percent_new_E = new_E / total_pop
    percent_new_I = new_I / total_pop
    percent_new_R = new_R / total_pop
    psum = percent_S + percent_E + percent_VD1 + percent_VD2 + percent_I + percent_R + percent_D

    assert 0.95 <= psum < 1.05, ("at time " + str(curr_time))

    return [int(sim_time), percent_S, percent_E, percent_VD1, percent_VD2, percent_I, percent_R, percent_new_E, percent_new_I, percent_new_R, percent_D, psum]

try:
    if __name__ == "__main__":
        # Read the data of all regions and their names
        with open(log_filename, "r") as log_file:
            # Read the file twice
            #   Once to get the total_pop
            #   A second time to calculate the data
            for i in range(2):
                log_file.seek(0, 0)
                line_num = 0

                # For each line, read a line then:
                for line in log_file:
                    # Strip leading and trailing spaces
                    line = line.strip()

                    # If a time marker is found that is not the current time
                    if line.isnumeric() and line != curr_time:
                        if curr_states and i == 1:
                            data.append(curr_states_to_df_row(curr_time, curr_states, sum(list(total_pop.values()))))

                        # Update new simulation time
                        curr_time = line
                        continue

                    # Create an re match objects from the current line
                    state_match = re.search(regex_state,line)
                    id_match    = re.search(regex_model_id,line)
                    if not (state_match and id_match):
                        continue

                    # Parse the state and id and insert into total_pop
                    cid     = id_match.group().lstrip('_')
                    state   = state_match.group().strip('<>')
                    state   = state.split(',')

                    if i == 0:
                        total_pop[cid] = float(state[0])
                    elif i == 1:
                        state = list(map(float, state))
                        curr_states[cid] = state

                    line_num += 1

            data.append(curr_states_to_df_row(curr_time, curr_states, sum(total_pop.values())))

        font = {"family" : "DejaVu Sans",
                "weight" : "normal",
                "size"   : 16}

        matplotlib.rc("font", **font)
        matplotlib.rc("lines", linewidth=2)

        try:
            os.makedirs(path)
        except OSError as error:
            raise error

        with open(base_name+"aggregate_timeseries.csv", 'w') as out_file:
            out_file.write("sim_time, S, E, VD1, VD2, I, R, New_E, New_I, New_R, D, pop_sum\n")
            for timestep in data:
                out_file.write(str(timestep).strip('[]')+"\n")

        columns = ["time", "susceptible", "exposed", "vaccinatedD1", "vaccinatedD2", "infected",
                    "recovered", "new_exposed", "new_infected", "new_recovered", "deaths", "error"]
        df_vis = pd.DataFrame(data, columns=columns)
        df_vis = df_vis.set_index("time")
        df_vis.to_csv(log_file_folder+"/states.csv")
        df_vis.head()
        x = list(df_vis.index)

        ### --- New EIR --- ###
        fig, ax = plt.subplots(figsize=(15,6))

        ax.plot(x, 100*df_vis["new_exposed"],   label="New exposed",   color=COLOR_EXPOSED,   linestyle=STYLE_EXPOSED)
        ax.plot(x, 100*df_vis["new_infected"],  label="New infected",  color=COLOR_INFECTED,  linestyle=STYLE_INFECTED)
        ax.plot(x, 100*df_vis["new_recovered"], label="New recovered", color=COLOR_RECOVERED, linestyle=STYLE_RECOVERED)
        plt.legend(loc="upper right")
        plt.title("Epidemic Aggregate New EIR Percentages")
        plt.xlabel("Time (days)")
        plt.ylabel("Population (%)")
        plt.savefig(base_name + "New_EIR.png")

        ### --- SEIRD/SEVIRD --- ###
        fig, axs = plt.subplots(2, figsize=(15,6))

        axs[0].plot(x, 100*df_vis["susceptible"], label="Susceptible", color=COLOR_SUSCEPTIBLE, linestyle=STYLE_SUSCEPTIBLE)
        if not (sum(df_vis['vaccinatedD1']) == 0 and sum(df_vis['vaccinatedD2']) == 0):
            axs[0].plot(x, 100*df_vis["vaccinatedD1"], label="Vaccinated D1", color=COLOR_DOSE1, linestyle=STYLE_DOSE1)
            axs[0].plot(x, 100*df_vis["vaccinatedD2"], label="Vaccinated D2", color=COLOR_DOSE2, linestyle=STYLE_DOSE2)
            axs[0].set_title("Epidemic Aggregate SEVIRD Percentages")
        else:
            axs[0].set_title("Epidemic Aggregate SEIR+D Percentages")

        axs[0].plot(x, 100*df_vis["exposed"],   label="Exposed",   color=COLOR_EXPOSED,   linestyle=STYLE_EXPOSED)
        axs[0].plot(x, 100*df_vis["infected"],  label="Infected",  color=COLOR_INFECTED,  linestyle=STYLE_INFECTED)
        axs[0].plot(x, 100*df_vis["recovered"], label="Recovered", color=COLOR_RECOVERED, linestyle=STYLE_RECOVERED)
        axs[0].set_ylabel("Population (%)")
        axs[0].legend(loc="upper right")

        axs[1].plot(x, 100*df_vis["deaths"],   label="Deceased", color=COLOR_DEAD,     linestyle=STYLE_DEAD)
        axs[1].plot(x, 100*df_vis["exposed"],  label="Exposed",  color=COLOR_EXPOSED,  linestyle=STYLE_EXPOSED)
        axs[1].plot(x, 100*df_vis["infected"], label="Infected", color=COLOR_INFECTED, linestyle=STYLE_INFECTED)
        axs[1].set_xlabel("Time (days)")
        axs[1].set_ylabel("Population (%)")
        axs[1].legend(loc="upper right")

        if not (sum(df_vis['vaccinatedD1']) == 0 and sum(df_vis['vaccinatedD2']) == 0):
            plt.savefig(base_name + "SEVIRD.png")
        else:
            plt.savefig(base_name + "SEIR+D.png")

        if not progress:
            print("\033[1;32mDone.\033[0m")
        else:
            done = True
            t.join()
    #if
except AssertionError as assertion:
    success = False
    done = True
    if progress:
        t.join()

    print("\n\033[31mASSERT:\033[0m 0.995 <= psum < 1.005", assertion)
    sys.exit(-1)
except KeyboardInterrupt as interrupt:
    success = False
    done = True
    if progress:
        t.join()

    print("\n\033[33mStopped by user\033[0m")
    sys.exit(-1)
except Exception as error:
    success = False
    done = True
    if progress:
        t.join()

    print("\n\033[31m" + str(error) + "\033[0m")
    sys.exit(-1)
