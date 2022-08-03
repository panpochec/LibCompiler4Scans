from pathlib import Path
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


def finding_rows(what, where, start, finish):
    """Function allowing scanning of a given matrix's in search for a given parameter.

    Parameters
    ----------
    where : list
        list to scan through
    what : str
        search for a given string
    start: int
        denominates begining of 'where'
    start: int
        denominates end of 'where'
    Returns
    -------
    list
        Returns list of where indexes numbers.
    """
    temp_list = []  # Temporary list, returned as output
    for p in range(0, len(where)):
        line_part = where[p]
        if str(what) in line_part[int(start): int(finish)]:
            temp_list.append(str(p))
        else:
            pass
    return temp_list


def finding_in_rows(what, which, where, start, finish, str_start, str_finish):
    """Function allowing scanning of a given matrix's in search for a given string.
    'str_start' and 'str_finish' allow searching only a given part of the line
    (useful when given letter or string is non-specific, like D denominating dihedrals
    but also denominating change of energy at the end of the line)

    Parameters
    ----------
    what : list
        search for a string from a given list
    which : int
        which of the 'what' elements to search for
    where : list
        Scanning of a given matrix
    start: int
        denominates begining of 'where'
    start: int
        denominates end of 'where'
    str_start : int
        allow searching only from a given point in a line of 'where'
    str_finish : int
        allow searching only to a given point in a line of 'where'

    Returns
    -------
    list
        Returns list of string contents.
    """
    temp_list = []  # Temporary list, returned as output
    for o in range(int(start), int(finish)):
        line_part = where[o]
        if str(what[which]) in line_part[str_start: str_finish]:
            temp_list.append(str(where[o]))
        else:
            pass
    return temp_list


def atomic_distance(where, atom1, atom2):
    """Function calculating distance between two atoms
        using xyz coordinates in a three-dimensional space

        Parameters
        ----------
        where : DataFrame
            xyz coordinates matrix (col1, col2, col3)
        atom1: int or string
            atom 1 number
        atom2: int or string
            atom 2 number

        Returns
        -------
        float
            Distance between atoms.
        """
    dist = math.sqrt(
                    (where.iloc[int(atom1), 0] - where.iloc[int(atom2), 0])
                    ** 2 +
                    (where.iloc[int(atom1), 1] - where.iloc[int(atom2), 1])
                    ** 2 +
                    (where.iloc[int(atom1), 2] - where.iloc[int(atom2), 2])
                    ** 2
                    )
    return dist

#TODO BridgeOne and BridgeTwo - done
#TODO BridgeOne and BridgeTwo energies - done
#TODO BridgeOne and BridgeTwo bridge atoms number - done
#TODO Ring atom numbers - done
#TODO Substitution pattern - done

#!!!!!!!!!!!!!!!!!!!!!! Building filelist !!!!!!!!!!!!!!!!!!!!!!

dirPath = Path.cwd()
result =  [f.name for f in dirPath.glob('*BridgeOne.log') if f.is_file()]
ile = len(result)

key = "SCF Done:"
key2 = "Optimized Parameters"
key3 = "Symbolic Z-matrix:"
key4 = "The following ModRedundant input section has been read:"
key5 = "Trust Radius"
key6 = " The following ModRedundant input section has been read:"
key7 = "Hirshfeld charges with hydrogens summed into heavy atoms"
key8 = "EQQ="
key9 = " Mulliken charges:"
key10 = "Sum of Mulliken charges ="
zrob = 0
df_list = []
df_parameters = []
bond = [0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4,
        0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85,
        0.9, 0.95, 1]
os.mkdir(r"./Graphs")

#!!!!!!!!!!!!!!!!!!!!!! Main file, iterating filelist !!!!!!!!!!!!!!!!!!!!!!

for y in range(0, ile):
    number = 0
    argu = result[y]
    argu2 = argu.replace("log", "csv")
    argu4 = argu.replace("_BridgeOne.log", "")
    argu5 = argu.replace("_BridgeOne.log", ".png")
    argu6 = argu.replace("BridgeOne", "BridgeTwo") # Accounting for files for both bridges
    argu7 = argu.replace("_BridgeOne", "")

    #!!!!!!!!!! BRIDGEONE !!!!!!!!!!!!

    # !!!!!!!!!! Searching for keyphrases associated with energies for
    file = open(argu, "r")
    lines = file.readlines()
    app = []
    app2 = []
    csv = []
    file.close()
    for number, line in enumerate(lines, 0):
        if key in line:
            app.append(number)
    for number, line in enumerate(lines, 0):
        if key2 in line:
            app2.append(number)
    arr = np.array(app)
    ile2 = len(app2)

    # Appending list with energies for each of the steps
    for x in range(0, ile2):
        num = int(app2[x])
        num2 = arr[arr < num].max()
        eneline = lines[num2]
        ene = eneline[26:40]
        ene2 = ene + " \n"
        csv.append(ene2)
#    LOGFILE = open(r"./CSV/" + argu2, "a")
#    LOGFILE.writelines(csv)

    # Acquiring bridge atoms numbers from modredundant
    for number, line in enumerate(lines, 0):
        if key6 in line:
            bridge_a_line = lines[number+1]
            break

    for x in range (0, 20):
        if bridge_a_line[x].isdigit():
            bridge_a_line = bridge_a_line.replace(bridge_a_line[0: bridge_a_line.index(bridge_a_line[x])], "")
            break
    bridge_a_oxygen = str(bridge_a_line[0:bridge_a_line.index(" ")])
    bridge_a_line = bridge_a_line.replace(bridge_a_oxygen, "")
    for x in range(0, 50):
        if bridge_a_line[0] == " ":
            bridge_a_line = bridge_a_line.replace(" ", "", 1)
        else:
            break
    bridge_a_hydrogen = str(bridge_a_line[0:bridge_a_line.index(" ")])
    bridge_a_line = bridge_a_line.replace(bridge_a_hydrogen, "")
    for x in range(0, 50):
        if bridge_a_line[0] == " ":
            bridge_a_line = bridge_a_line.replace(" ", "", 1)
        else:
            break
    bridge_a_nitrogene = str(bridge_a_line[0:bridge_a_line.index(" ")])

    # !!!!!!!!!!!!! BRIDGETWO !!!!!!!!!!!!!

    # !!!!!!!!! Searching for keyphrases associated with energies
    file = open(argu6, "r")
    lines_two = file.readlines()
    app = []
    app2 = []
    csv2 = []
    file.close()
    for number, line in enumerate(lines_two, 0):
        if key in line:
            app.append(number)
    for number, line in enumerate(lines_two, 0):
        if key2 in line:
            app2.append(number)
    arr = np.array(app)
    ile2 = len(app2)

    # !!!!!!!!!!! Appending list with energies for each of the steps
    for x in range(0, ile2):
        num = int(app2[x])
        num2 = arr[arr < num].max()
        eneline = lines_two[num2]
        ene = eneline[26:40]
        ene2 = ene + " \n"
        csv2.append(ene2)
    #    LOGFILE = open(r"./CSV/" + argu2, "a")
    #    LOGFILE.writelines(csv)

    # Acquiring bridge atoms numbers from modredundant
    for number, line in enumerate(lines, 0):
        if key6 in line:
            bridge_b_line = lines_two[number + 1]
            break

    for x in range(0, 20):
        if bridge_b_line[x].isdigit():
            bridge_b_line = bridge_b_line.replace(bridge_b_line[0: bridge_b_line.index(bridge_b_line[x])], "")
            break
    bridge_b_oxygen = str(bridge_b_line[0:bridge_b_line.index(" ")])
    bridge_b_line = bridge_b_line.replace(bridge_b_oxygen, "")
    for x in range(0, 50):
        if bridge_b_line[0] == " ":
            bridge_b_line = bridge_b_line.replace(" ", "", 1)
        else:
            break
    bridge_b_hydrogen = str(bridge_b_line[0:bridge_b_line.index(" ")])
    bridge_b_line = bridge_b_line.replace(bridge_b_hydrogen, "")
    for x in range(0, 50):
        if bridge_b_line[0] == " ":
            bridge_b_line = bridge_b_line.replace(" ", "", 1)
        else:
            break
    bridge_b_nitrogene = str(bridge_b_line[0:bridge_b_line.index(" ")])

    # !!!!!! FULL LIST of atoms and parameters !!!!!!

    # Parameter matrix and atom-type matrix
    for number, line in enumerate(lines, 0):  # Line before atom-type matrix
        if key3 in line:
            a_matrix_start = int(number) + 1
    for number, line in enumerate(lines, 0):  # Line after atom-type matrix
        if key4 in line:
            a_matrix_end = int(number) - 1
            p_matrix_start = int(number) + 13  # Line begining parameter matrix
            break
    atom_matrix = lines[a_matrix_start:a_matrix_end] # Creating atoms matrix
    atom_matrix[0] = "X"   ### to take
    print("Atoms matrix")
    print(atom_matrix)

    for number, line in enumerate(lines, 0):  # Finding line after parameter matrix
        if key5 in line:
            p_matrix_end = int(number) -1
            break
    parameter_matrix = lines[p_matrix_start:p_matrix_end]  # Parameter Matrix
    print("Parameter matrix:")
    print(parameter_matrix)

    # Finding dihedral rows  in parameter matrix
    dihedral_matrix = finding_rows("D", parameter_matrix, 0, 20)

    # Finding indexes for rows containing distances in par. matrix
    distances_row_list = finding_rows("R", parameter_matrix, 0, len(parameter_matrix))

    # Finding indexes for rows containing distances in par. matrix
    angles_row_list = finding_rows("A", parameter_matrix, 0, len(parameter_matrix))

    # Finding oxygens in rows
    rows_for_O = finding_rows("O", atom_matrix, 0, -1)
    for x in range(0, len(rows_for_O)):
        rows_for_O[x] = "s" + rows_for_O[x] + "s"
    # Finding nitrogens in rows
    rows_for_N = finding_rows("N", atom_matrix, 0, -1)
    for x in range(0, len(rows_for_N)):
        rows_for_N[x] = "s" + rows_for_N[x] + "s"
    # Finding hydrogens in rows
    rows_for_H = finding_rows("H", atom_matrix, 0, -1)
    for x in range(0, len(rows_for_H)):
        rows_for_H[x] = "s" + rows_for_H[x] + "s"
    # Finding carbons in rows
    rows_for_C = finding_rows("C", atom_matrix, 0, -1)
    for x in range(0, len(rows_for_C)):
        rows_for_C[x] = "s" + rows_for_C[x] + "s"
    # Finding bromines in rows
    rows_for_Br = finding_rows("Br", atom_matrix, 0, -1)
    for x in range(0, len(rows_for_Br)):
        rows_for_Br[x] = "s" + rows_for_Br[x] + "s"

    non_bridge_nitrogenes = []
    for x in range(0, len(rows_for_N)):
        if "s" + bridge_a_nitrogene + "s" not in rows_for_N[x]\
                and "s" + bridge_b_nitrogene + "s" not in rows_for_N[x]:
            non_bridge_nitrogenes.append(rows_for_N[x])

    # Hirshfeld charges list
    file = open(r"./CHARGES/" + argu7, "r")
    lines_three = file.readlines()
    file.close()
    for number, line in enumerate(lines_three, 0):
        if key7 in line:
            start_charges = number + 2
    for number, line in enumerate(lines_three, 0):
        if key8 in line:
            end_charges = number - 1

    for number, line in enumerate(lines_three, 0):
        if key9 in line:
            start_M_charges = number + 2
    for number, line in enumerate(lines_three, 0):
        if key10 in line:
            end_M_charges = number - 1

    charges_list = lines_three[start_charges: end_charges]
    for x in range(0, len(charges_list)):
        charges_list[x] = str(charges_list[x])[12:]
        if str(charges_list[x])[0] == " ":
            charges_list[x] = str(charges_list[x])[1:]
    for x in range(0, len(charges_list)):
        charges_list[x] = str(charges_list[x])[0:str(charges_list[x]).index(" ")]
    t = 0
    charges_list_clean = []
    for x in range(1, len(atom_matrix)):
        if "s" + str(x) + "s" in rows_for_H:
            charges_list_clean.append("Hydrogen")
        else:
            charges_list_clean.append(charges_list[t])
            t += 1

    charges_M_list = lines_three[start_M_charges: end_M_charges]
    for x in range(0, len(charges_M_list)):
        charges_M_list[x] = float(charges_M_list[x][12:])

    # DataFrame for dihedral angles

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    col5 = []

    # Due to the construction of lines in parameter matrix,
    # searching for atom numbers and building dataframe is done
    # using a trick:
    # We iterate thorugh each line of data matrix searching for a specific
    # denominator, then appending respective column list and deleting the appended
    # part with the denominator to make further iterations easier.
    # In the case of dihedrals :
    # ! D44   D(16,17,39,40)        179.9269         estimate D2E/DX2
    # First we delete the part "! D44   D("
    # Then we take count to the first ",", and append column list with 16
    # and use the newly created entry '16' and denominator ',' and
    # delete them from line. The resulting line:
    # 17,39,40)        179.9269         estimate D2E/DX2
    # We repeat the process until we acquire all of the parameters.
    # Same process is used for angles and distances

    for x in range(int(dihedral_matrix[0]), int(dihedral_matrix[-1])+1):
        help_string = str(parameter_matrix[x])
        help_string = help_string.replace(str(parameter_matrix[x])[0:11], "")
        col1.append(int(str(help_string)[0: str(help_string).index(",", 0, 50)]))
        help_string = help_string.replace(str(col1[x - (int(dihedral_matrix[0]))]) + ",", "", 1)
        col2.append(int(str(help_string)[0: str(help_string).index(",", 0)]))
        help_string = help_string.replace(str(col2[x - (int(dihedral_matrix[0]))]) + ",", "", 1)
        col3.append(int(str(help_string)[0: str(help_string).index(",", 0)]))
        help_string = help_string.replace(str(col3[x - (int(dihedral_matrix[0]))]) + ",", "", 1)
        col4.append(int(str(help_string)[0: str(help_string).index(")", 0)]))
        help_string = help_string.replace(str(col4[x - (int(dihedral_matrix[0]))]) + ")", "", 1)
        help_string = help_string.replace(" ", "")
        col5.append(float(str(help_string)[0: str(help_string).index("e", 0)]))

    df_dihedral = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Atom 3': col3,
        'Atom 4': col4,
        'Value': col5
    })

    col1 = []
    col2 = []
    col3 = []
    col4 = []
    # DataFrame for angles
    for x in range(int(angles_row_list[0]), int(angles_row_list[-1])+1):
        help_string = str(parameter_matrix[x])
        help_string = help_string.replace(str(parameter_matrix[x])[0:11], "")
        col1.append(int(str(help_string)[0: str(help_string).index(",", 0, 50)]))
        help_string = help_string.replace(str(col1[x - (int(angles_row_list[0]))]) + ",", "", 1)
        col2.append(int(str(help_string)[0: str(help_string).index(",", 0)]))
        help_string = help_string.replace(str(col2[x - (int(angles_row_list[0]))]) + ",", "", 1)
        col3.append(int(str(help_string)[0: str(help_string).index(")", 0)]))
        help_string = help_string.replace(str(col3[x - (int(angles_row_list[0]))]) + ")", "", 1)
        help_string = help_string.replace(" ", "")
        if "Fro" in help_string:
            col4.append(float(str(help_string)[0: str(help_string).index("F", 0)]))
        else:
            col4.append(float(str(help_string)[0: str(help_string).index("e", 0)]))


    df_angle = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Atom 3': col3,
        'Value': col4
    })

    # Building DataFrame for distances
    # Auxilary data
    col1 =[]
    col2 = []
    col3 = []
    for x in range(int(distances_row_list[0]), int(distances_row_list[-1])+1):
        col1.append(int(str(parameter_matrix[x])[11: int(str(parameter_matrix[x]).index(','))]))
        col2.append(int(str(parameter_matrix[x])[int(str(
            parameter_matrix[x]).index(',')) + 1: int(str(parameter_matrix[x]).index(')'))]))
        col3.append(float(str(parameter_matrix[x])[33: int(str(parameter_matrix[x]).index(' ', 35, 45))]))

    # DataFrame
    distance_matrix = pd.DataFrame({
        'Atom 1': col1,
        'Atom 2': col2,
        'Atom 3': col3
    })

    # Finding bridge carbons (carbon in each ring closest to bridge oxygen)
    l = -1
    for row in distance_matrix.iterrows():
        l += 1
        if " " + bridge_a_oxygen + ".000" in str(row):
            if str(distance_matrix.iloc[l, 0]) != bridge_a_hydrogen \
                    and str(distance_matrix.iloc[l, 0]) != bridge_a_oxygen:
                bridge_a_carbon = str(distance_matrix.iloc[l, 0])
            elif str(distance_matrix.iloc[l, 1]) != bridge_a_hydrogen \
                    and str(distance_matrix.iloc[l, 1]) != bridge_a_oxygen:
                bridge_a_carbon = str(distance_matrix.iloc[l, 1])

    x = -1
    for rowz in distance_matrix.iterrows():
        x += 1
        if " " + bridge_b_oxygen + ".00" in str(rowz):
            if str(distance_matrix.iloc[x, 0]) != bridge_b_hydrogen \
                    and str(distance_matrix.iloc[x, 0]) != bridge_b_oxygen:
                bridge_b_carbon = str(distance_matrix.iloc[x, 0])
            elif str(distance_matrix.iloc[x, 1]) != bridge_b_hydrogen \
                    and str(distance_matrix.iloc[x, 1]) != bridge_b_oxygen:
                bridge_b_carbon = str(distance_matrix.iloc[x, 1])

    # Dihedral values
    dihedral_value_list = []
    for x in range(int(dihedral_matrix[0]), int(dihedral_matrix[-1])):
        b = str(parameter_matrix[x])
        dihedral_value_list.append(b[29:41])

    for x in range(0, len(dihedral_value_list)):
        dihedral_value_list[x] = float(str(dihedral_value_list[x]).replace(" ", ""))

    dihedral_value_diff = []  # Dihedral values correction
    for x in range(0, len(dihedral_value_list)):
        if abs(dihedral_value_list[x]) < 90:
            dihedral_value_diff.append(abs(dihedral_value_list[x]))
        else:
            dihedral_value_diff.append(abs(dihedral_value_list[x]) - 180)

    # Molecular Mass Calculator

    mass_count = ''.join(atom_matrix)  # Joining list to make it countable
    carbons = int(mass_count.count("C"))
    oxygens = int(mass_count.count("O"))
    bromines = int(mass_count.count("B"))
    hydrogens = int(mass_count.count("H"))
    nitrogenes = int(mass_count.count("N"))
    final_mass = carbons*12 + oxygens*16 + bromines*80 + hydrogens + nitrogenes*14
    print("Molecular mass")
    print(final_mass)


    # Conversion of energies to kcal for DataFrames

    # BridgeOne

    # Building new data from list csv
    # to make them floats, as well
    # as diff. from 1st point and
    # conversion to kcal/mol
    csv3 = []
    for kk in range(0, len(csv)):
        csv3.append((float(csv[kk].replace("\n", ""))
                     - float(csv[0].replace("\n", "")))*627.51)

    # BridgeTwo #
    csv4 = []
    for kk in range(0, len(csv2)):
        csv4.append((float(csv2[kk].replace("\n", ""))
                     - float(csv2[0].replace("\n", "")))*627.51)

    # Compound group assessment
    if "meta" in argu:
        compound_type = "meta"
    elif "para" in argu:
        compound_type = "para"
    elif "orto" in argu:
        compound_type = "orto"
    else:
        compound_type = "Not specified"

    if "pcm" in argu4:
        model = "Polarizable Continuum Model"
        model_operator = 'pcm'
    else:
        model = "Gas Phase"
        model_operator = ''

    # Buliding naming scheme for ring atoms

    # Bridge One

    # Building carbon list for ring A using
    # dihedral angles with bridge carbon and eleminating any containing
    # anything besides carbons, then appending list with
    # entries from the search.
    carbons_ring_a_helper = []
    z = -1
    for row1 in df_dihedral.iterrows(): # Iterating rows
        z += 1
        if " " + bridge_a_carbon + ".0" in str(row1):
            if "s" + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_H and "s" \
                    + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_Br and "s" \
                    + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_N and "s" \
                    + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_O:
                if "s" + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_H and "s" \
                        + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_Br and "s" \
                        + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_N and "s" \
                        + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_O:
                    if "s" + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_H and "s" \
                            + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_Br and "s" \
                            + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_N and "s" \
                            + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_O:
                        if "s" + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_H and "s" \
                                + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_Br and "s" \
                                + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_N and "s" \
                                + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_O:
                            carbons_ring_a_helper.append(df_dihedral.iloc[z, 0])
                            carbons_ring_a_helper.append(df_dihedral.iloc[z, 1])
                            carbons_ring_a_helper.append(df_dihedral.iloc[z, 2])
                            carbons_ring_a_helper.append(df_dihedral.iloc[z, 3])

    # Eliminating multiple entries
    carbons_ring_a_temp = []
    for x in range(0, len(carbons_ring_a_helper)):
        if carbons_ring_a_helper[x] not in carbons_ring_a_temp:
            carbons_ring_a_temp.append(carbons_ring_a_helper[x])

    # Eliminating non-bridge carbon using common distance with bridge nitrogen
    for x in range(0, len(carbons_ring_a_temp)):
        for row in distance_matrix.iterrows():
            if " " + str(carbons_ring_a_temp[x]) + ".00" in str(row) and " " \
                    + str(bridge_a_nitrogene) + ".00" in str(row):
                non_ring_carbon_a = carbons_ring_a_temp[x]
    carbons_ring_a = [x for x in carbons_ring_a_temp if x is not non_ring_carbon_a]

    # Assigning specific positional value

    carbon_ring_a_5 = int(bridge_a_carbon)
    carbons_ring_a.remove(carbon_ring_a_5)

    # Searching list of ring carbons for one with common distance with
    # bridge nitrogene, then eliminating from entry from list
    for x in range(0, len(carbons_ring_a)):
        for tup in df_angle.iterrows():
            if " " + str(carbons_ring_a[x]) + ".00" in str(tup) and " " \
                    + str(bridge_a_nitrogene) + ".00" in str(tup):
                carbon_ring_a_6 = carbons_ring_a[x]
    carbons_ring_a.remove(carbon_ring_a_6)

    # Searching list of ring carbons for one with common dihedral with
    # bridge nitrogene and without carbon 5, then eliminating entry from list
    for x in range(0, len(carbons_ring_a)):
        for row in df_dihedral.iterrows():
            if " " + str(carbons_ring_a[x]) + ".00" in str(row) and " " \
                    + str(bridge_a_nitrogene) + ".00" in str(row) and " " \
                    + str(carbon_ring_a_5) + ".00" not in str(row):
                carbon_ring_a_1 = carbons_ring_a[x]
    carbons_ring_a.remove(carbon_ring_a_1)

    # For the rest, searching list of ring carbons for one with common dihedral with
    # other carbons but always excluding some, for unambiguity then eliminating entry from list
    for x in range(0, len(carbons_ring_a)):
        for row in df_dihedral.iterrows():
            if " " + str(carbons_ring_a[x]) + ".00" in str(row) and " " \
                    + str(carbon_ring_a_1) + ".00" in str(row) and " " \
                    + str(carbon_ring_a_6) + ".00" in str(row) and " " \
                    + str(non_ring_carbon_a) + ".00" in str(row) and " " \
                    + str(carbon_ring_a_5) + ".00" not in str(row):
                carbon_ring_a_2 = carbons_ring_a[x]
    carbons_ring_a.remove(carbon_ring_a_2)
    for x in range(0, len(carbons_ring_a)):
        for row in df_dihedral.iterrows():
            if " " + str(carbons_ring_a[x]) + ".00" in str(row) and " " \
                    + str(carbon_ring_a_1) + ".00" in str(row) and " " \
                    + str(carbon_ring_a_6) + ".00" in str(row) and " " + str(carbon_ring_a_2) \
                    + ".00" in str(row) and " " + str(carbon_ring_a_5) \
                    + ".00" not in str(row):
                carbon_ring_a_3 = carbons_ring_a[x]
    carbons_ring_a.remove(carbon_ring_a_3)
    for row in df_dihedral.iterrows():
        if " " + str(carbons_ring_a[0]) + ".00" in str(row) and " " \
                + str(carbon_ring_a_1) + ".00" in str(row) and " " \
                + str(carbon_ring_a_6) + ".00" not in str(row) and " " \
                + str(carbon_ring_a_2) + ".00" in str(row) and " " \
                + str(carbon_ring_a_3) + ".00" in str(row):
            carbon_ring_a_4 = carbons_ring_a[0]

    ring_a_substitution_sites = [carbon_ring_a_1, carbon_ring_a_2, carbon_ring_a_3, carbon_ring_a_4]

    # Bridge Two ###

    # Building carbon list for ring A using
    # dihedral angles with bridge carbon and eleminating any containing
    # anything besides carbons, then appending list with
    # entries from the search.
    carbons_ring_b_helper = []
    z = -1
    for row2 in df_dihedral.iterrows():  # Iterating rows
        z += 1
        if " " + bridge_b_carbon + ".00" in str(row2):
            if "s" + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_H and "s" \
                    + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_Br and "s" \
                    + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_N and "s" \
                    + str(df_dihedral.iloc[z, 0]) + "s" not in rows_for_O:
                if "s" + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_H and "s" \
                        + str(df_dihedral.iloc[z, 1]) + "s" not in rows_for_Br and "s" \
                        + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_N and "s" \
                        + str(df_dihedral.iloc[z, 1]) + "s"not in rows_for_O:
                    if "s" + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_H and "s" \
                            + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_Br and "s" \
                            + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_N and "s" \
                            + str(df_dihedral.iloc[z, 2]) + "s" not in rows_for_O:
                        if "s" + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_H and "s" \
                                + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_Br and "s" \
                                + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_N and "s" \
                                + str(df_dihedral.iloc[z, 3]) + "s" not in rows_for_O:
                            carbons_ring_b_helper.append(df_dihedral.iloc[z, 0])
                            carbons_ring_b_helper.append(df_dihedral.iloc[z, 1])
                            carbons_ring_b_helper.append(df_dihedral.iloc[z, 2])
                            carbons_ring_b_helper.append(df_dihedral.iloc[z, 3])

    # Eliminating multiple entries
    carbons_ring_b_temp = []
    for x in range(0, len(carbons_ring_b_helper)):
        if carbons_ring_b_helper[x] not in carbons_ring_b_temp:
            carbons_ring_b_temp.append(carbons_ring_b_helper[x])

    # Eliminating non-bridge carbon using common distance with bridge nitrogen
    for x in range(0, len(carbons_ring_b_temp)):
        for row in distance_matrix.iterrows():
            if " " + str(carbons_ring_b_temp[x]) + ".00" in str(row) and " " \
                    + str(bridge_b_nitrogene) + ".00" in str(row):
                non_ring_carbon_b = carbons_ring_b_temp[x]
    carbons_ring_b = [x for x in carbons_ring_b_temp if x is not non_ring_carbon_b]

    # Assigning specific positional value

    carbon_ring_b_5 = int(bridge_b_carbon)
    carbons_ring_b.remove(carbon_ring_b_5)

    # Searching list of ring carbons for one with common distance with
    # bridge nitrogene, then eliminating from entry from list
    for x in range(0, len(carbons_ring_b)):
        for tup in df_angle.iterrows():
            if " " + str(carbons_ring_b[x]) + ".00" in str(tup) and " " \
                    + str(bridge_b_nitrogene) + ".00" in str(tup):
                carbon_ring_b_6 = carbons_ring_b[x]
    carbons_ring_b.remove(carbon_ring_b_6)

    # Searching list of ring carbons for one with common dihedral with
    # bridge nitrogene and without carbon 5, then eliminating entry from list
    for x in range(0, len(carbons_ring_b)):
        for row in df_dihedral.iterrows():
            if " " + str(carbons_ring_b[x]) + ".00" in str(row) and " " \
                    + str(bridge_b_nitrogene) + ".00" in str(row) and " " \
                    + str(carbon_ring_b_5) + ".00" not in str(row):
                carbon_ring_b_1 = carbons_ring_b[x]
    carbons_ring_b.remove(carbon_ring_b_1)

    # For the rest, searching list of ring carbons for one with common dihedral with
    # other carbons but always excluding some, for unambiguity then eliminating entry from list
    for x in range(0, len(carbons_ring_b)):
        for row in df_dihedral.iterrows():
            if " " + str(carbons_ring_b[x]) + ".00" in str(row) and " " \
                    + str(carbon_ring_b_1) + ".00" in str(
                    row) and " " + str(carbon_ring_b_6) + ".00" in str(row) and " " + str(non_ring_carbon_b) \
                    + ".00" in str(
                    row) and " " + str(carbon_ring_b_5) + ".00" not in str(row):
                carbon_ring_b_2 = carbons_ring_b[x]
    carbons_ring_b.remove(carbon_ring_b_2)
    for x in range(0, len(carbons_ring_b)):
        for row in df_dihedral.iterrows():
            if " " + str(carbons_ring_b[x]) + ".00" in str(row) and " " \
                    + str(carbon_ring_b_1) + ".00" in str(
                    row) and " " + str(carbon_ring_b_6) + ".00" in str(row) and " " + str(carbon_ring_b_2) + ".00" in str(
                    row) and " " + str(carbon_ring_b_5) + ".00" not in str(row):
                carbon_ring_b_3 = carbons_ring_b[x]
    carbons_ring_b.remove(carbon_ring_b_3)
    for row in df_dihedral.iterrows():
        if " " + str(carbons_ring_b[0]) + ".00" in str(row) and " " \
                + str(carbon_ring_b_1) + ".00" in str(
                row) and " " + str(carbon_ring_b_6) + ".00" not in str(row) and " " + str(carbon_ring_b_2) + ".00" in str(
                row) and " " + str(carbon_ring_b_3) + ".00" in str(row):
            carbon_ring_b_4 = carbons_ring_b[0]

    ring_b_substitution_sites = [carbon_ring_b_1, carbon_ring_b_2, carbon_ring_b_3, carbon_ring_b_4]

    # Connector atoms in middle rings

    for row in distance_matrix.iterrows():
        for x in range(0, len(rows_for_C)):
            if " " + str(bridge_b_nitrogene) + ".00" in str(row) and " " \
                    + str(rows_for_C[x].replace("s", "")) + ".00" in str(row) and " " \
                    + str(non_ring_carbon_b) + ".00" not in str(row):
                ring_b_connector = rows_for_C[x].replace("s", "")
            elif " " + str(bridge_a_nitrogene) + ".00" in str(row) and " " \
                    + str(rows_for_C[x].replace("s", "")) + ".00" in str(row) and " " \
                    + str(non_ring_carbon_a) + ".00" not in str(row):
                ring_a_connector = rows_for_C[x].replace("s", "")

    # Specifing substitution pattern
    substitution_name = argu4.replace(compound_type, "")
    substitution_name = substitution_name.replace("pcm", "")
    if substitution_name == "":
        substitution_group = ["Not applicable", "Not applicable"]
        ring_substitution_postion = ["Not applicable", "Not applicable"]
        which_ring = ["Not applicable", "Not applicable"]
        substitution_pattern = ["Not applicable", "Not applicable"]
        a_group = "None"
        a_pattern = "None"
        b_group = "None"
        b_pattern = "None"
        index_name = compound_type + model_operator
        a_substitution_id = "None,None"
        b_substitution_id = "None,None"

    else:
        for x in range(0, len(substitution_name)):
            if substitution_name[x].isdigit():
                substitution_group_crude = substitution_name[0:x]
                substitution_pattern_crude = substitution_name[x:]
                break
        substitution_group = []
        ring_substitution_postion = []
        if len(substitution_group_crude) == 4:
            substitution_group.append(substitution_group_crude[0:2])
            substitution_group.append(substitution_group_crude[2:])
        else:
            substitution_group.append(substitution_group_crude[0:2])

        if len(rows_for_Br) + len(non_bridge_nitrogenes) == 2:
            substitution_type = "di"
        else:
            ring_substitution_postion.append(substitution_pattern_crude[0])
            substitution_type = "mono"

        if substitution_type == "di" and len(substitution_group_crude) == 2:
            substitution_group.append(substitution_group[0])

        for x in range(0, len(substitution_group)):
            if substitution_group[x] == "br":
                substitution_group[x] = "Br"
            elif substitution_group[x] == "nh":
                substitution_group[x] = 'NH2'
            elif substitution_group[x] == "no":
                substitution_group[x] = 'NO2'

        # Checking if non-bridge nitrogen is in amine or nitro group
        amine_nitrogene = []
        nitro_nitrogene = []
        nitro_nitrogene_help = []
        for x in range(0, len(non_bridge_nitrogenes)):
            for row in distance_matrix.iterrows():
                if " " + str(non_bridge_nitrogenes[x].replace("s", "")) + ".000" in str(row):
                    for z in range(0, len(rows_for_O)):
                        if " " + str(rows_for_O[z].replace("s", "")) + ".000" in str(row):
                            nitro_nitrogene_help.append(non_bridge_nitrogenes[x])

        for x in range(0, len(nitro_nitrogene_help)):
            if nitro_nitrogene_help[x] not in nitro_nitrogene:
                nitro_nitrogene.append(nitro_nitrogene_help[x])

        for x in range(0, len(non_bridge_nitrogenes)):
            if non_bridge_nitrogenes[x] in nitro_nitrogene:
                pass
            else:
                amine_nitrogene.append(non_bridge_nitrogenes[x])

        substitution_pattern = []
        which_ring = []
        b = 0
        nh = 0
        no = 0
        # Assigining substituent to substitution site
        for x in range(0, len(substitution_group)):
            if substitution_group[x] == "Br":
                for row in distance_matrix.iterrows():
                    for c in range(0, len(ring_a_substitution_sites)):
                        if " " + rows_for_Br[b].replace("s", "") \
                                + ".0" in str(row) and " " \
                                + str(ring_a_substitution_sites[c]) + ".0" in str(row):
                           substitution_pattern.append(str(c+1))
                           which_ring.append("A")
                        elif " " + rows_for_Br[b].replace("s", "") + ".0" in str(row) and " " \
                                + str(ring_b_substitution_sites[c]) \
                                + ".0" in str(row):
                            substitution_pattern.append(str(c + 1))
                            which_ring.append("B")
                b += 1
            elif substitution_group[x] == "NH2":
                for row in distance_matrix.iterrows():
                    for c in range(0, len(ring_a_substitution_sites)):
                        if " " + amine_nitrogene[nh].replace("s", "") \
                                + ".0" in str(row) and " " + str(ring_a_substitution_sites[c]) + ".0" in str(row):
                           substitution_pattern.append(str(c+1))
                           which_ring.append("A")
                        elif " " + amine_nitrogene[nh].replace("s", "") \
                                + ".0" in str(row) and " " + str(ring_b_substitution_sites[c]) \
                                + ".0" in str(row):
                            substitution_pattern.append(str(c + 1))
                            which_ring.append("B")
                nh += 1
            elif substitution_group[x] == "NO2":
                for row in distance_matrix.iterrows():
                    for c in range(0, len(ring_a_substitution_sites)):
                        if " " + nitro_nitrogene[no].replace("s", "") \
                                + ".00" in str(row) and " " + str(ring_a_substitution_sites[c]) \
                                + ".00" in str(row):
                           substitution_pattern.append(str(c+1))
                           which_ring.append("A")
                        elif " " + nitro_nitrogene[no].replace("s", "") \
                                + ".00" in str(row) and " " + str(ring_b_substitution_sites[c]) \
                                + ".00" in str(row):
                            substitution_pattern.append(str(c + 1))
                            which_ring.append("B")
                no += 1
        a_group = "None"
        a_pattern = "None"
        b_group = "None"
        b_pattern = "None"

        for x in range(0, len(which_ring)):
            if which_ring[x] == "A":
                a_group = substitution_group[x]
                a_pattern = substitution_pattern[x]
            elif which_ring[x] == "B":
                b_group = substitution_group[x]
                b_pattern = substitution_pattern[x]

        # Writing index for final dataframe
        if len(substitution_group) == 2:
            index_name = compound_type + substitution_group[0] + \
                         str(substitution_pattern[0]) + substitution_group[1] + \
                         str(substitution_pattern[1]) + model_operator
        elif len(substitution_group) == 1:
            index_name = compound_type + substitution_group[0] + str(substitution_pattern[0]) + model_operator

        a_substitution_id = a_group + "," + a_pattern
        b_substitution_id = b_group + "," + b_pattern

    x = 0
    for row in distance_matrix.iterrows():
        if " " + str(bridge_a_oxygen) + ".00" in str(row) and " " + str(bridge_a_hydrogen) + ".00" in str(row):
            bridge_a_length = distance_matrix.iloc[x, 2]
        x += 1

    x = 0
    for row in distance_matrix.iterrows():
        if " " + str(bridge_b_oxygen) + ".00" in str(row) and " " + str(bridge_b_hydrogen) + ".00" in str(row):
            bridge_b_length = distance_matrix.iloc[x, 2]
        x += 1


    # Finding activation energy and second minimum
    for max in range(0, len(csv3)):
        if csv3[max] > csv3[max+1]:
            activation_energy_bridgeone = csv3[max]
            break
        else:
            pass

    for min in range(csv3.index(activation_energy_bridgeone), len(csv3)):
        if csv3[min] < csv3[min+1]:
            second_minimum_bridgeone = csv3[min]
            break
        else:
            pass
    for max in range(0, len(csv4)):
        if csv4[max] > csv4[max+1]:
            activation_energy_bridgetwo = csv4[max]
            break
        else:
            pass

    for min1 in range(csv4.index(activation_energy_bridgetwo), len(csv4)):
        if csv4[min1] < csv4[min1+1]:
            second_minimum_bridgetwo = csv4[min1]
            break
        else:
            pass

    index_atom = ['x']
    colx = [0.0]
    coly = [0.0]
    colz = [0.0]
    for x in range(1, len(atom_matrix)):
        index_atom.append(str(atom_matrix[x])[0:2])
        colx.append(float(str(atom_matrix[x])[22:30]))
        coly.append(float(str(atom_matrix[x])[32:41]))
        colz.append(float(str(atom_matrix[x])[42:51]))

    atom_postions = pd.DataFrame({
        'X': colx,
        'Y': coly,
        'Z': colz
    }, index=index_atom)

    bridge_a_NOdist = atomic_distance(atom_postions, bridge_a_nitrogene, bridge_a_oxygen)
    bridge_a_NHdist = atomic_distance(atom_postions, bridge_a_nitrogene, bridge_a_hydrogen)
    bridge_a_CNdist = atomic_distance(atom_postions, bridge_a_nitrogene, carbon_ring_a_1)
    bridge_a_COdist = atomic_distance(atom_postions, carbon_ring_a_4, bridge_a_oxygen)
    bridge_b_NOdist = atomic_distance(atom_postions, bridge_b_nitrogene, bridge_b_oxygen)
    bridge_b_NHdist = atomic_distance(atom_postions, bridge_b_nitrogene, bridge_b_hydrogen)
    bridge_b_CNdist = atomic_distance(atom_postions, bridge_b_nitrogene, carbon_ring_b_1)
    bridge_b_COdist = atomic_distance(atom_postions, carbon_ring_b_4, bridge_b_oxygen)

    # Building Data Frame for a single compound
    df_part = pd.DataFrame({
        'bond length diff.': bond[0: len(csv3)],
        'energy': csv3,
        'header_name': argu,
        'header_mass': final_mass,
        'header_type': compound_type,
        'header_activation_ene': activation_energy_bridgeone,
        'header_second_min': second_minimum_bridgeone
    })
    df_list.append(df_part)
    df_part2 = pd.DataFrame({
        'bond length diff.': bond[0: len(csv4)],
        'energy': csv4,
        'header_name': argu,
        'header_mass': final_mass,
        'header_type': compound_type,
        'header_activation_ene': activation_energy_bridgetwo,
        'header_second_min': second_minimum_bridgetwo
    })
    df_list.append(df_part2)
    # Appending longfile for R use

    df_part_parameters = pd.DataFrame({
        'File name': argu4,
        'Mass': final_mass,
        'Type': compound_type,
        'Model': model,
        'Ring A substituent': a_group,
        'Ring A substituted in postion': a_pattern,
        'Ring A starting lenght': bridge_a_length,
        'Activation energy - bridge A': activation_energy_bridgeone,
        'Second minimum - bridge A': second_minimum_bridgeone,
        'Ring B substituent': b_group,
        'Ring B substituted in postion': b_pattern,
        'Ring B starting lenght': bridge_b_length,
        'Activation energy - bridge B': activation_energy_bridgetwo,
        'Second minimum - bridge B': second_minimum_bridgetwo,
        'A ring carbon 5 Hirschfeld charge': charges_list_clean[int(carbon_ring_a_5)-1],
        'A ring carbon 6 Hirschfeld charge': charges_list_clean[int(carbon_ring_a_6) - 1],
        'A ring carbon 4 Hirschfeld charge': charges_list_clean[int(carbon_ring_a_4) - 1],
        'A ring carbon 3 Hirschfeld charge': charges_list_clean[int(carbon_ring_a_3) - 1],
        'A ring carbon 2 Hirschfeld charge': charges_list_clean[int(carbon_ring_a_2) - 1],
        'A ring carbon 1 Hirschfeld charge': charges_list_clean[int(carbon_ring_a_1) - 1],
        'A non-ring carbon Hirschfeld charge': charges_list_clean[int(non_ring_carbon_a) - 1],
        'Connector carbon A Hirschfeld charge': charges_list_clean[int(ring_a_connector)-1],
        'B ring carbon 5 Hirschfeld charge': charges_list_clean[int(carbon_ring_b_5)-1],
        'B ring carbon 6 Hirschfeld charge': charges_list_clean[int(carbon_ring_b_6) - 1],
        'B ring carbon 4 Hirschfeld charge': charges_list_clean[int(carbon_ring_b_4) - 1],
        'B ring carbon 3 Hirschfeld charge': charges_list_clean[int(carbon_ring_b_3) - 1],
        'B ring carbon 2 Hirschfeld charge': charges_list_clean[int(carbon_ring_b_2) - 1],
        'B ring carbon 1 Hirschfeld charge': charges_list_clean[int(carbon_ring_b_1) - 1],
        'B non-ring carbon Hirschfeld charge': charges_list_clean[int(non_ring_carbon_b) - 1],
        'Connector carbon B Hirschfeld charge': charges_list_clean[int(ring_b_connector)-1],
        'Bridge A nitrogen charge': charges_list_clean[int(bridge_a_nitrogene) - 1],
        'Bridge A oxygen charge': charges_list_clean[int(bridge_a_oxygen) - 1],
        'Bridge B nitrogen charge': charges_list_clean[int(bridge_b_nitrogene) - 1],
        'Bridge B oxygen charge': charges_list_clean[int(bridge_b_oxygen) - 1],
        'A ring carbon 5 Mulliken charge': charges_M_list[int(carbon_ring_a_5) - 1],
        'A ring carbon 6 Mulliken charge': charges_M_list[int(carbon_ring_a_6) - 1],
        'A ring carbon 4 Mulliken charge': charges_M_list[int(carbon_ring_a_4) - 1],
        'A ring carbon 3 Mulliken charge': charges_M_list[int(carbon_ring_a_3) - 1],
        'A ring carbon 2 Mulliken charge': charges_M_list[int(carbon_ring_a_2) - 1],
        'A ring carbon 1 Mulliken charge': charges_M_list[int(carbon_ring_a_1) - 1],
        'A non-ring carbon Mulliken charge': charges_M_list[int(non_ring_carbon_a) - 1],
        'Connector carbon A Mulliken charge': charges_M_list[int(ring_a_connector) - 1],
        'B ring carbon 5 Mulliken charge': charges_M_list[int(carbon_ring_b_5) - 1],
        'B ring carbon 6 Mulliken charge': charges_M_list[int(carbon_ring_b_6) - 1],
        'B ring carbon 4 Mulliken charge': charges_M_list[int(carbon_ring_b_4) - 1],
        'B ring carbon 3 Mulliken charge': charges_M_list[int(carbon_ring_b_3) - 1],
        'B ring carbon 2 Mulliken charge': charges_M_list[int(carbon_ring_b_2) - 1],
        'B ring carbon 1 Mulliken charge': charges_M_list[int(carbon_ring_b_1) - 1],
        'B non-ring carbon Mulliken charge': charges_M_list[int(non_ring_carbon_b) - 1],
        'Connector carbon B Mulliken charge': charges_M_list[int(ring_b_connector) - 1],
        'Bridge A nitrogen Mulliken charge': charges_M_list[int(bridge_a_nitrogene) - 1],
        'Bridge A oxygen Mulliken charge': charges_M_list[int(bridge_a_oxygen) - 1],
        'Bridge B nitrogen Mulliken charge': charges_M_list[int(bridge_b_nitrogene) - 1],
        'Bridge B oxygen Mulliken charge': charges_M_list[int(bridge_b_oxygen) - 1],
        'Bridge A nitrogen - oxygen distance': bridge_a_NOdist,
        'Bridge A nitrogen - hydrogen distance': bridge_a_NHdist,
        'Bridge A nitrogen - carbon 1 distance': bridge_a_CNdist,
        'Bridge A carbon 4  - oxygen distance': bridge_a_COdist,
        'Bridge B nitrogen - oxygen distance': bridge_b_NOdist,
        'Bridge B nitrogen - hydrogen distance': bridge_b_NHdist,
        'Bridge B nitrogen - carbon 1 distance': bridge_b_CNdist,
        'Bridge B carbon 4  - oxygen distance': bridge_b_COdist,
        'A ring substitution ID': a_substitution_id,
        'B ring substitution ID': b_substitution_id

    }, index=[index_name])

    #df_part.plot(x='bond length diff.', y='energy', kind='scatter', title=argu4 ,
    #xlabel =r'Bond length diff. [$\AA$]', ylabel =r'$\delta$E[kcal/mol]', ylim=[-0.5,16.5] )

    # Ploting for a single compound
    fig, axs = plt.subplots(2, 1, figsize=(5, 7))
    fig.suptitle("Difference in energy while scanning ")
    axs[0].set_title('Bridge A, substitution: ' + a_pattern + ", " + a_group)
    axs[0].plot(df_part.iloc[:, 0], df_part.iloc[:, 1], 'bo')
    axs[0].set(ylabel=r'$\delta$E [kcal/mol]', xlim=[-0.05, 0.8], ylim=[-0.5, 11])
    axs[0].grid(axis="y")
    axs[1].set_title('Bridge B, substitution: ' + b_pattern + ", " + b_group)
    axs[1].plot(df_part2.iloc[:, 0], df_part2.iloc[:, 1], 'ro')
    axs[1].set(xlabel=r'Bond length diff. [$\AA$]', ylabel=r'$\delta$E [kcal/mol]', xlim=[-0.05, 0.8], ylim=[-0.5, 11])
    axs[1].grid(axis="y")
    plt.subplots_adjust(left=0.2,
                        bottom=0.1,
                        right=0.9,
                        top=0.9,
                        wspace=0.3,
                        hspace=0.4)
    #plt.show()
    df_parameters.append(df_part_parameters)
    plt.savefig(r"./Graphs/" + argu5)
pd.concat(df_parameters).to_csv("ScanLibrary.csv", index=True)         ###Joining final DataFrame to .csv
print("FileDone")
