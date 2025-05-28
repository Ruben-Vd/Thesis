import tkinter as tk
import numpy as np
import math
from tkinter import scrolledtext

##Fixed parameters
#These parameters are fixed for the calculation, they are not given by the user
fric_coeff = 0.15 #friction coefficient.
pressure_angle = 20 #in degrees.
Y_DT = 1.0 #set to 1 as it takes into account tolerance classes and that is out of the scope of this thesis.
K = 1.35 #The total product of the K factors are set to 1.35.

# a bunch of lists, which are used to store all the different values
Epsilon_alpha_2 = [] #transverse contact ratio in stage 2
Epsilon_alpha_1 = [] #transverse contact ratio in stage 1
Epsilon_alpha__1 = [] #transverse contact ratio in stage 1'
Epsilon_alpha_0 = [] #transverse contact ratio in stage 0
B2  = [] #width of stage 2
B1  = [] #width of stage 1
B_1  = [] #width of stage 1'
B0  = [] #width of stage 0
HA_1 = []  #helix angle in stage 1
HA_2 =[] #helix angle in stage 2
HA_0 = [] # helix angle in stage 0
HA__1 = []  #helix angle in stage 1'
A = [] #center distances in basic wolfrom so [[center distance in stage 1 for first solution, center distance in stage 2 for first solution], [center distance in stage 1 for second solution, center distance in stage 2 for second solution], ...]
A_pre = [] #same as A but for the pregearing
f_R2 =[]  # loss factor in stage 2 BW
f_R1 = []  # loss factor in stage 1 BW
f_S = []   # loss factor in stage 0 pre
f_P0 = []   # loss factor in stage 1' pre
Y_ST = 2 #set to 2 as in the norm
Y_X = 1 #as in norm
Z_W = 1 #as in norm
Z_X = 1 #as in norm


clearance = 1 #mm !!!! used in the neighbouring condition, 1mm clearance on the diameter of the gear, so 0.5mm on each side
clearance_PS = 1.5 #mm used in the determination of the maximal sum of profile shifts possible, actually another slightly severe neighbouring condition
margin = 0.5  #used in the extra neighboring condition                                                                                                                                                      
Material_properties = {"Steel 20MNCr-5": [206000000000, 0.3, 1200, 850, 425,  0.003, 0.957, 1500, 4.8], "Steel 18CrNiMo7-6": [206000000000, 0.3, 1200, 850, 500, 0.003, 0.957, 1500, 4.8],"Steel 16 MnCrS5":[206000000, 0.3, 1000, 695, 425, 0.003, 0.957, 4.8], "Steel 40CrMoV13-9":[206000000, 0.3, 950, 750, 420, 0.1005, 0.99, 1250, 16]} #E-modulus, poisson ratio, tensile strength in MPA, yield point in MPA, sigma_F_limiet, rho', YRrelT, sigma_H_limiet, Rz (mean to valley roughness)  kijk naar KISSsoft materialen

Tooth_thickness_allowances = {"DIN 3967 cd25": [-0.054, -0.084, -0.07, -0.11], "DIN 3967 cd26": [-0.054, -0.1040, -0.07, -0.13],"DIN 3967 f22":[-0.014, -0.022, -0.019, -0.029]} # tooth thickness allowances for the different types of gears, cd25, cd26 and f22 are the different norms of DIN, the first number is the upper tolerance for the external gear, the second the lower for the external gear, the third the upper for the internal gear and the fourth teh lower for the internal gear 


output = tk.Tk() # Create the main window
output.title("Gearbox Calculation Results") # Set the title of the window
output.geometry("900x550+1000+150")  # Initial size of the window 900 is width, 550 is height, +1000+150 is the position of the window on the screen

 
def clear_output(): # function that will clear the output text area when the button is pressed
    output_text.delete("1.0", tk.END)

def clear_input(): # function that will clear the input fields when the button is pressed
    for entry in entry_boxes: # Clear all entry fields
        entry.delete(0, tk.END)
    
    allowance_var.set("DIN 3967 cd25")  # Reset allowance dropdown
    material_var.set("Steel 20MNCr-5")  # Reset material dropdown
    gear_var.set('Spur')
    hole_var.set("No")  # Reset hole dropdown

    # Hide elements that were dynamically shown
    config_label.grid_remove()
    config_entry.grid_remove()
    config_label_pre.grid_remove()
    config_entry_pre.grid_remove()
    button4.grid_remove()
    button3.grid(row=10, column=4, columnspan=2, padx=10, pady=10, sticky="s") #put calculate button again in the right bottom corner

    # Also hide the hole diameter input if it was shown
    hole_diameter_label.grid_remove()
    hole_diameter_entry.grid_remove()

    Epsilon_alpha_2 = [] #transverse contact ratio in stage 2
    Epsilon_alpha_1 = [] #transverse contact ratio in stage 1
    Epsilon_alpha__1 = [] #transverse contact ratio in stage 1'
    Epsilon_alpha_0 = [] #transverse contact ratio in stage 0
    B2  = [] #width of stage 2
    B1  = [] #width of stage 1
    B_1  = [] #width of stage 1'
    B0  = [] #width of stage 0
    HA_1 = []  #helix angle in stage 1
    HA_2 =[] #helix angle in stage 2
    HA_0 = [] # helix angle in stage 0
    HA__1 = []  #helix angle in stage 1'
    A = [] #center distances in basic wolfrom so [[center distance in stage 1 for first solution, center distance in stage 2 for first solution], [center distance in stage 1 for second solution, center distance in stage 2 for second solution], ...]
    A_pre = [] #same as A but for the pregearing
    f_R2 =[]  # loss factor in stage 2 BW
    f_R1 = []  # loss factor in stage 1 BW
    f_S = []   # loss factor in stage 0 pre
    f_P0 = []   # loss factor in stage 1' pre

# Add a button above the text area to clear the output
button_frame = tk.Frame(output)
button_frame.grid(row=0, column=0, padx=10, pady=5, sticky="new")  # "new" to keep it at the top

# Add "Clear Output" button inside the button frame
clear_button = tk.Button(button_frame, text="Clear Output", command=clear_output)
clear_input_button = tk.Button(button_frame, text="Reset input screen", command=clear_input)

clear_button.grid(row=0, column=0, padx=5, pady=5)
clear_input_button.grid(row=0, column=1, padx=5, pady=5)


# Create a frame to hold the text and scrollbars
frame = tk.Frame(output)
frame.grid(row=1, column=0, padx=10, pady=10, sticky="nsew")

# Configure the main window to resize both vertically and horizontally
output.grid_rowconfigure(0, weight=0)
output.grid_rowconfigure(1, weight=1)  # Row 0 (frame) should expand vertically
output.grid_columnconfigure(0, weight=1)  # Column 0 (frame) should expand horizontally

# Create the ScrolledText widget
output_text = scrolledtext.ScrolledText(frame, width=108, height=30, wrap=tk.NONE)
output_text.grid(row=0, column=0, sticky="ew")  # Use sticky to ensure resizing in both directions

# Add a horizontal scrollbar below the text area
h_scrollbar = tk.Scrollbar(frame, orient=tk.HORIZONTAL, command=output_text.xview)
h_scrollbar.grid(row=1, column=0, sticky="ew")  # Make sure the horizontal scrollbar expands 

# Configure the ScrolledText widget to work with the horizontal scrollbar
output_text.config(xscrollcommand=h_scrollbar.set)
output_text.tag_configure("bold", font=("Helvetica", 10, "bold"))

# Make sure the scrollbar doesn't affect the layout too much
frame.grid_rowconfigure(0, weight=1)  # The first row (ScrolledText) should expand (weight = 1)
frame.grid_columnconfigure(0, weight=1)  # The second row (horizontal scrollbar) should also expand

# Allow the entire window to resize in both directions
output_text.grid(sticky="nsew")

#message to the user when starting the program
output_text.insert(tk.END, "In this window, the results of the program will be shown.\n", 'bold')
output_text.insert(tk.END, "\n")
output_text.insert(tk.END, "Please fill in all the input fields and select the desired type of gears, material and cable hole diameter (if any).\n", 'bold')
output_text.insert(tk.END, "\n")
output_text.insert(tk.END, "When one wants to use helical gears, the transverse module must be given!\n", 'bold')
output_text.insert(tk.END, "\n")
output_text.insert(tk.END, "Once you filled in all the values, click on the 'Preliminary study' button in the bottom left corner to start.\n", 'bold')

def SF_min(m_n): #minimum safety factor against bending according to ISO 6336 (KISSsoft)
    if m_n <= 0.5:
        SF_min = 0.6
    elif m_n==1.0:
        SF_min = 1.2
    elif m_n >=2.0:
        SF_min = 1.4
    elif 0.5 < m_n <1.0: #linear interpolating
        SF_min = 0.6 + (m_n - 0.5)*1.2
    elif 1 < m_n < 2:
        SF_min =  1.2 + (m_n - 1.0)*0.2
    return SF_min

def SH_min(m_n): #minimum safety factor against contact stress according to ISO 6336 (KISSsoft)
    if m_n <= 0.5:
        SH_min = 0.6
    elif m_n == 1.0:
        SH_min = 0.9
    elif m_n >= 2.0:
        SH_min = 1.0
    elif 0.5 < m_n < 1.0: #linear interpolating
        SH_min = 0.6 + (m_n - 0.5) * (0.9 - 0.6) / (1.0 - 0.5)
    elif 1.0 < m_n < 2.0:
        SH_min = 0.9 + (m_n - 1.0) * (1.0 - 0.9) / (2.0 - 1.0)
    return SH_min

def inverse_involute(theta, tol=1e-8, max_iter=20): #to find the value of alpha_wt from an involute value
    phi = np.sqrt(2 * theta)  # Initial guess
    for _ in range(max_iter):
        f = np.tan(phi) - phi - theta
        df = 1 / np.cos(phi)**2 - 1
        phi_new = phi - f / df
        if abs(phi_new - phi) < tol:
            return phi_new
        phi = phi_new
    raise ValueError("Newton-Raphson did not converge")

def calculate_preliminary(): #this function tries to find a minimum module that is needed in order to be strong enough in theory
    global m_min_bending_stage_1, m_min_bending_stage_2, module1_ok, module2_ok, module0ok, module_1ok, m_min_bending_stage_0, m_min_bending_stage__1 #will be used in the calculate preliminary and calculate basic wolfrom function to use the correct module
    output_text.delete("1.0", tk.END) #delete text in output field to replace it with the outcome of this function
    ## Get the values from the entry boxes and convert them to the correct types
    module1_ok = True #from the start is it assumed to have the right (high enough) values 
    module2_ok = True
    module_1ok = True
    module0ok = True
    values = [entry.get() for entry in entry_boxes]
    max_radius = (float(values[0])/2)*10**(-3)
    T_out = float(values[1])
    axial_length_wolf = float(values[3])
    axial_length_pre = float(values[10])
    iw = float(values[4])
    i0 = float(values[11])
    NOP_pre = int(values[13])
    NOP_wolf = int(values[6])
    eta_w = int(values[9])/100
    eta_0 = int(values[16])/100
    m_1_orig = float(values[7]) #helical module for helical gear
    m_2_orig = float(values[8])
    m__1_orig = float(values[15])
    m_0_orig =float(values[14])
    max_radius_sun = max_radius + 3.5*10*10**(-3) + 1.25*10*10**(-3) - 0.5*10**(-3) - 11*0.5*10**(-3)
    material = material_var.get()
    max_bending = Material_properties[material][3]
    type_of_gear = gear_var.get()
    tangential_forces = []
    axial_forces = []
    radial_forces = []
    normal_forces = []

    # Calculate the tangential forces for the different gears
    F_T_R2 = T_out / max_radius
    F_T_P2 = F_T_R2 / NOP_wolf
    F_T_R1 = (-(1 - 1 / (iw * eta_w)) * T_out) / max_radius
    F_T_P1 = F_T_R1 / NOP_wolf
    F_T_s = -(1 / (iw * i0)) * (1 / (eta_w * eta_0)) * T_out / max_radius_sun
    F_T_P0 = F_T_s / NOP_pre
    F_T_R_1 = ((1 / (iw * i0)) * (1 / (eta_w * eta_0)) + (1 / (iw * eta_w))) * T_out / max_radius
    F_T_P_1 = F_T_R_1 / NOP_pre

    # Add to list in order
    forces = [F_T_R2, F_T_P2, F_T_R1, F_T_P1, F_T_s, F_T_P0, F_T_R_1, F_T_P_1]
    for F in forces:
        tangential_forces.append(abs(round(F, 2)))

    #Y factors for internal and external gears
    Y_B_int = 1
    Y_B_ext = 1

    Y_F_ext = 2.37 #https://sdp-si.com/resources/elements-of-metric-gear-technology/page8.php#Section17, 
    Y_S_ext = 1 #set to one in order to certainly not overestimate the module

    Y_F_int = 2.06
    Y_S_int = 1 #set to one in order to certainly not overestimate the module

    m_1 = m_1_orig
    m_2 = m_2_orig
    m_0 = m_0_orig
    m__1 = m__1_orig
    if type_of_gear == 'Spur': #depending on the type of gears, the calculation is different, as there is a dependence on the module in the width calculation of the helical gears
        b1,b2 = calculate_width(F_T_P1, F_T_P2, axial_length_wolf, type_of_gear, m_1, m_2)
        b0,b_1 = calculate_width(F_T_P0, F_T_P_1, axial_length_pre, type_of_gear, m_0, m__1)
        helix_angle_input_wolf = 0
        helix_angle_output_wolf = 0
        Y_beta_output_wolf = 1
        Y_beta_input_wolf =  1
        helix_angle_input_pre = 0
        helix_angle_output_pre = 0
        Y_beta_output_pre = 1
        Y_beta_input_pre =  1
        m_min_bending_stage_2 = min( abs(F_T_P2)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_output_wolf*K/(b2*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_wolf)),  abs(F_T_R2) * Y_F_int*Y_DT*Y_S_int*Y_B_int*Y_beta_output_wolf*K/(b2*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_wolf)))
        m_min_bending_stage_1 = min( abs(F_T_P1)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_input_wolf*K/(b1*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_wolf)), abs(F_T_R1) * Y_F_int*Y_DT*Y_S_int*Y_B_int*Y_beta_input_wolf*K/(b1*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_wolf)))
        m_min_bending_stage_0 = min( abs(F_T_P0)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_input_pre*K/(b0*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_pre)),  abs(F_T_s) * Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_input_pre*K/(b0*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_pre)))
        m_min_bending_stage__1 = min( abs(F_T_P_1)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_output_pre*K/(b_1*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_pre)), abs(F_T_R_1) * Y_F_int*Y_DT*Y_S_int*Y_B_int*Y_beta_output_pre*K/(b_1*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_pre)))
    else:
        if type_of_gear == "Helical, eps_beta = 1":
            overlap_ratio = 1
        if type_of_gear == "Helical, eps_beta = 2":
            overlap_ratio = 2
        error__1 = 5
        error_0 = 5
        error_1 = 5
        error_2 = 5
        count = 0
        while error_1 and error_2 and error_0 and error__1> 10**(-1): #try to converge to a certain minimum module becuase here it can not be avoided to use the inserted moduli
            if count > 500: #when the number of iterations is too high, no convergence is reached and the loop is stopped
                m_min_bending_stage_0 = np.nan
                m_min_bending_stage_1 = np.nan
                m_min_bending_stage__1 = np.nan
                m_min_bending_stage_2 = np.nan
                break
            b1,b2 = calculate_width(F_T_P1, F_T_P2, axial_length_wolf, type_of_gear, m_1, m_2)
            b0,b_1 = calculate_width(F_T_P0, F_T_P_1, axial_length_pre, type_of_gear, m_0, m__1)
            helix_angle_input_wolf = np.arctan(overlap_ratio*np.pi*m_1/b1)
            helix_angle_output_wolf = np.arctan(overlap_ratio*np.pi*m_2/b2)
            Y_beta_input_wolf =  (1-overlap_ratio*helix_angle_input_wolf*180/(np.pi*120))*1/((np.cos(helix_angle_input_wolf))**3)
            Y_beta_output_wolf =  (1-overlap_ratio*helix_angle_output_wolf*180/(np.pi*120))*1/((np.cos(helix_angle_output_wolf))**3)
            helix_angle_input_pre = np.arctan(overlap_ratio*np.pi*m_0/b0)
            helix_angle_output_pre = np.arctan(overlap_ratio*np.pi*m__1/b_1)
            Y_beta_input_pre =  (1-overlap_ratio*helix_angle_input_pre*180/(np.pi*120))*1/((np.cos(helix_angle_input_pre))**3)
            Y_beta_output_pre =  (1-overlap_ratio*helix_angle_output_pre*180/(np.pi*120))*1/((np.cos(helix_angle_output_pre))**3)
            m_min_bending_stage_2 = min( abs(F_T_P2)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_output_wolf*K/(b2*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_wolf)),  abs(F_T_R2) * Y_F_int*Y_DT*Y_S_int*Y_B_int*Y_beta_output_wolf*K/(b2*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_wolf)))
            m_min_bending_stage_1 = min( abs(F_T_P1)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_input_wolf*K/(b1*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_wolf)), abs(F_T_R1) * Y_F_int*Y_DT*Y_S_int*Y_B_int*Y_beta_input_wolf*K/(b1*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_wolf)))
            m_min_bending_stage_0 = min( abs(F_T_P0)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_input_pre*K/(b0*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_pre)),  abs(F_T_s) * Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_input_pre*K/(b0*10**(-3)*max_bending*10**6*np.cos(helix_angle_input_pre)))
            m_min_bending_stage__1 = min( abs(F_T_P_1)*Y_F_ext*Y_DT*Y_S_ext*Y_B_ext*Y_beta_output_pre*K/(b_1*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_pre)), abs(F_T_R_1) * Y_F_int*Y_DT*Y_S_int*Y_B_int*Y_beta_output_pre*K/(b_1*10**(-3)*max_bending*10**6*np.cos(helix_angle_output_pre)))
            error_0 = abs(m_0-m_min_bending_stage_0*10**3)
            error__1 = abs(m__1-m_min_bending_stage__1*10**3)
            error_1 = abs(m_1-m_min_bending_stage_1*10**3)
            error_2 = abs(m_2-m_min_bending_stage_2*10**3)
            m_1 = m_min_bending_stage_1*10**3
            m_2 = m_min_bending_stage_2*10**3
            m_0 = m_min_bending_stage_0*10**3
            m__1 = m_min_bending_stage__1*10**3
            count +=1
        m_1 = m_1_orig #finally set the values of the tangential moduli back to the desired ones in order to compare.
        m_2 = m_2_orig
        m_0 = m_0_orig
        m__1 = m__1_orig

    for i in range(8):
        # Select the correct helix angle
        if i in [0, 1]:
            helix_angle = helix_angle_output_wolf * 180 / np.pi
        elif i in [2, 3]:
            helix_angle = helix_angle_input_wolf * 180 / np.pi
        elif i in [4, 5]:
            helix_angle = helix_angle_input_pre * 180 / np.pi
        elif i in [6, 7]:
            helix_angle = helix_angle_output_pre * 180 / np.pi

        # Append calculated forces
        axial_forces.append(calculate_F_a(tangential_forces[i], helix_angle))
        radial_forces.append(calculate_F_r(tangential_forces[i], pressure_angle, helix_angle))
        normal_forces.append(calculate_F_n(tangential_forces[i], pressure_angle))

    output_text.insert(tk.END, "Results of the preliminary study:\n", "bold")
    output_text.insert(tk.END, " ", "bold")
    # Formatting output into a table-like structure
    output_text.insert(tk.END, "{:<25}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format(" ", "R2", "P2", "R1", "P1", "s", "P0", "R1'", "P1'")) 
    output_text.insert(tk.END, "{:<25}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_t_min [N]", round(F_T_R2, 2), round(F_T_P2, 2), round(F_T_R1, 2), round(F_T_P1, 2), round(F_T_s, 2), round(F_T_P0, 2), round(F_T_R_1, 2), round(F_T_P_1, 2)))
    output_text.insert(tk.END, "{:<25}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_a_min [N]", round(axial_forces[0], 2), round(axial_forces[1], 2), round(axial_forces[2], 2), round(axial_forces[3], 2), round(axial_forces[4], 2), round(axial_forces[5], 2), round(axial_forces[6], 2),round(axial_forces[7], 2)))
    output_text.insert(tk.END, "{:<25}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_r_min [N]", round(radial_forces[0], 2), round(radial_forces[1], 2), round(radial_forces[2], 2), round(radial_forces[3], 2), round(radial_forces[4], 2), round(radial_forces[5], 2), round(radial_forces[6], 2), round(radial_forces[7], 2)))
    output_text.insert(tk.END, "{:<25}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_n_min [N]", round(normal_forces[0], 2), round(normal_forces[1], 2), round(normal_forces[2], 2), round(normal_forces[3], 2), round(normal_forces[4], 2), round(normal_forces[5], 2), round(normal_forces[6], 2), round(normal_forces[7], 2)))

    # Add calculated values
    output_text.insert(tk.END, "\nMinimum module required:\n")
    output_text.insert(tk.END, "Minimum (tangential) module required in input stage of basic Wolfrom: {:.6f} mm\n".format(m_min_bending_stage_1 * 10**3))
    output_text.insert(tk.END, "Minimum (tangential) module required in output stage of basic Wolfom: {:.6f} mm\n".format(m_min_bending_stage_2 * 10**3))
    output_text.insert(tk.END, "Minimum (tangential) module required in input stage of pregearing: {:.6f} mm\n".format(m_min_bending_stage_0 * 10**3))
    output_text.insert(tk.END, "Minimum (tangential) module required in output stage of pregearing: {:.6f} mm\n".format(m_min_bending_stage__1 * 10**3))

    # Check if the modules are acceptable
    output_text.tag_configure("bold_red", font=("Helvetica", 10, "bold"), foreground="red")

    # Check if the modules are acceptable and print warnings in bold
    if np.isnan(m_min_bending_stage_0):
        output_text.insert(tk.END, "Preliminary study did not converge for the input stage of the pregearing, the given module will further be used. It is up to the user to check if this module is too low or not.\n", "bold_red")
    if np.isnan(m_min_bending_stage_1):
        output_text.insert(tk.END, "Preliminary study did not converge for the input stage of the basic Wolfrom, the given module will further be used. It is up to the user to check if this module is too low or not.\n", "bold_red")
    if np.isnan(m_min_bending_stage__1):
        output_text.insert(tk.END, "Preliminary study did not converge for the output stage of the pregearing, the given module will further be used. It is up to the user to check if this module is too low or not.\n", "bold_red")
    if np.isnan(m_min_bending_stage_2):
        output_text.insert(tk.END, "Preliminary study did not converge for the output stage of the basic Wolfrom, the given module will further be used. It is up to the user to check if this module is too low or not.\n", "bold_red")    
    
    if m_1 < m_min_bending_stage_1 * 10**3:
        module1_ok = False #if the flag is set to false, the module will be adapted in the calculate_basic_wolfrom function
        output_text.insert(tk.END, "Warning: Input stage module in basic Wolfrom is too low!\n", "bold_red")

    if m_2 < m_min_bending_stage_2 * 10**3:
        module2_ok = False
        output_text.insert(tk.END, "Warning: Output stage module in basic Wolfrom is too low!\n", "bold_red")

    if m_0 < m_min_bending_stage_0 * 10**3:
        module0ok = False #if the flag is set to false, the module will be adapted in the calculate_pregearing function
        output_text.insert(tk.END, "Warning: Input stage module of pregearing is too low!\n", "bold_red")

    if m__1 < m_min_bending_stage__1 * 10**3:
        module_1ok = False
        output_text.insert(tk.END, "Warning: Output stage module of pregearing is too low!\n", "bold_red")
    output_text.insert(tk.END, "\n")
    output_text.insert(tk.END, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", "bold")
    output_text.insert(tk.END, "Click on the middle button to start the calculation for the basic Wolfrom part.", "bold")
    output.mainloop()    

def calculate_Y_beta_ISO(overlap_ratio, helix_angle): # to calculate helix angle factor, the input of np.cos must be in radians, calculated according to the ISO 6336:3 norm
    if overlap_ratio > 1:
        overlap_ratio = 1
    if helix_angle>30:
        helix_angle = 30
    return (1-overlap_ratio*helix_angle/120)*1/np.cos(helix_angle*np.pi/180)**3 

def calculate_Y_F_ISO(Kind_of_gear, NOT, helix_angle, prof_shift_coeff_generated,m_n, alpha_n, rho_fP, Epsilon_alpha_, x): 
    #according to ISO 6336:3 2019, profile shift coefficient is with the generating PS taken into account, x is the 'real' profile shift
    if Kind_of_gear == 'ext':
        T = np.pi/3
        d = NOT*m_n/np.cos(helix_angle*np.pi/180) 
        hfp = 1.25*m_n
        beta_b = np.arcsin(np.sin(helix_angle*np.pi/180)*np.cos(alpha_n*np.pi/180))
        eps_alpha_n = Epsilon_alpha_/(np.cos(beta_b)**2)
        if overlap_ratio == 0 and eps_alpha_n < 2:
            f_eps = 1
        if overlap_ratio == 0 and eps_alpha_n >= 2:
            f_eps = 0.7
        if 0<overlap_ratio<1 and eps_alpha_n < 2:
            f_eps = (1 - overlap_ratio + overlap_ratio/eps_alpha_n)**(1/2)
        if 0 < overlap_ratio < 1 and eps_alpha_n >=2:
            f_eps = ((1-overlap_ratio)/2 + overlap_ratio/eps_alpha_n)**(1/2)
        if overlap_ratio>=1:
            f_eps = eps_alpha_n**(-0.5)
        Zn = NOT/(np.cos(beta_b)**2 * np.cos(helix_angle*np.pi/180))
        d_n = m_n*Zn
        dbn = d_n*np.cos(alpha_n*np.pi/180)
        da = d + 2*m_n*(1+x)
        dan = d_n + da - d
        den = 2* NOT/abs(NOT)*np.sqrt((np.sqrt(( (dan/2)**2 - (dbn/2)**2) ) - (np.pi*d*np.cos(helix_angle*np.pi/180)*np.cos(alpha_n*np.pi/180))/(abs(NOT))*(eps_alpha_n-1))**2 + (dbn/2)**2)   
        alpha_en = np.arccos(dbn/den)
        gamma_e = (0.5*np.pi + 2 * prof_shift_coeff_generated* np.tan(alpha_n*np.pi/180))/Zn + np.tan(alpha_n*np.pi/180)-alpha_n*np.pi/180 - np.tan(alpha_en) + alpha_en
        alpha_Fen = alpha_en - gamma_e
        E = np.pi/4*m_n - hfp*np.tan(alpha_n*np.pi/180) - (1-np.sin(alpha_n*np.pi/180))*rho_fP/np.cos(alpha_n*np.pi/180)  
        G = rho_fP/m_n - hfp/m_n + prof_shift_coeff_generated
        H = 2/Zn*(np.pi/2 - E/m_n) - T
        iteration = 0
        theta_init = np.pi/6
        while iteration <5:
            theta = 2*G/Zn*np.tan(theta_init)-H
            theta_init = theta
            iteration +=1
        sFn = m_n*Zn*np.sin(np.pi/3-theta)+ m_n*np.sqrt(3)*(G/np.cos(theta)-rho_fP/m_n)
        rho_F = m_n*(rho_fP/m_n + 2*G**2/(np.cos(theta)*(Zn*np.cos(theta)**2-2*G)))
        hFe = m_n/2*((np.cos(gamma_e) - np.sin(gamma_e)*np.tan(alpha_Fen))*den/m_n - Zn*np.cos(np.pi/3 - theta) - (G/np.cos(theta) - rho_fP/m_n))
        Y_F = (6*hFe/m_n * np.cos(alpha_Fen))/((sFn/m_n)**2 * np.cos(alpha_n*np.pi/180))*f_eps
    else:
        NOT = -NOT
        if abs(NOT) > 60: #the same code as KISSsoft uses, found by having contact with them.
            z0 = int(min(50, abs(NOT) * 0.33))
        elif abs(NOT) > 11:
            if abs(NOT) >19:
                z0 = int(max(12, abs(NOT) * 0.66))
            elif abs(NOT) > 16:
                z0 = int(max(10, abs(NOT)*0.66))
            else:
                z0 = abs(NOT) - 2
        else:
            z0 = abs(NOT)
        theta = np.pi/3
        x0 = 0
        haP0 = 1.25*m_n
        rho_a0 = 0.38*m_n
        beta_b = math.asin(math.sin(helix_angle*math.pi/180) * math.cos(alpha_n*math.pi/180))
        Z0v = z0/(math.cos(helix_angle*math.pi/180) * math.cos(beta_b)**2)
        Zn = NOT/(math.cos(beta_b)**2*math.cos(helix_angle*math.pi/180))
        eps_alphan = Epsilon_alpha_/(math.cos(beta_b)**2)
        ksi = 2*(x0+prof_shift_coeff_generated)/(Z0v+Zn)*math.tan(alpha_n*math.pi/180) + math.tan(alpha_n*math.pi/180) - alpha_n*math.pi/180
        alpha_w0_init = np.cbrt(3*ksi)
        counter = 0
        while counter <3:
            alpha_w0 = alpha_w0_init + (ksi - math.tan(alpha_w0_init)+alpha_w0_init)/(math.tan(alpha_w0_init)**2)
            alpha_w0_init = alpha_w0
            counter +=1
        alpha_w0 = alpha_w0_init
        a0 = m_n*(Z0v+Zn)/2*(math.cos(alpha_n*math.pi/180))/(math.cos(alpha_w0))
        u0 = Z0v/Zn
        rw = a0/(1+u0)
        rw0 = rw*u0
        rb0 = 0.5*m_n*Z0v*math.cos(alpha_n*math.pi/180)
        rM = m_n*(Z0v/2 + haP0/m_n +x0 - rho_a0/m_n)
        alpha_M = math.acos(rb0/rM)
        delta_alpha = (0.5*math.pi+2*x0*math.tan(alpha_n*math.pi/180))/Z0v -rho_a0/rb0 + math.tan(alpha_n*math.pi/180) - alpha_n*math.pi/180 - math.tan(alpha_M) + alpha_M
        YM = rM*math.cos(delta_alpha)
        psi_0 = math.pi/Zn +theta
        psi_init = psi_0
        lamba = rw0/rM *math.cos(psi_init)
        y = psi_init - math.acos(lamba) + delta_alpha + (psi_init-psi_0)/u0
        y_ = 1+1/u0 -rw0/rM*(math.sin(psi_init)/(math.sqrt(1-lamba**2)))
        while abs(y/y_) > 10**(-6):
            psi = psi_init -y/y_
            lamba = rw0/rM *math.cos(psi)
            y = psi - math.acos(lamba) + delta_alpha + (psi-psi_0)/u0
            y_ = 1+1/u0 -rw0/rM*(math.sin(psi)/(math.sqrt(1-lamba**2)))
            psi_init =psi
        psi = psi_init
        delta = math.acos(rw0/rM*math.cos(psi))
        omega_0 = delta - psi - delta_alpha
        deltah_ = YM - rw0*math.cos(omega_0)
        delta_h = (deltah_ *math.sin(psi))/(math.sin(psi+omega_0))
        K = delta_h/math.sin(psi)
        X = rw*math.sin(psi-theta) - (K+rho_a0)*math.cos(theta)
        Y = rw*math.cos(psi-theta) - (K+rho_a0)*math.sin(theta)
        sFn = 2*X
        rho_F = K**2/((rw0*rw)/(rw0+rw) * math.sin(psi)+K) + rho_a0
        d = NOT * m_n/np.cos(helix_angle*np.pi/180)
        dn = d/(math.cos(beta_b)**2)
        dbn = dn*math.cos(alpha_n*math.pi/180)
        ha = m_n*(1+x)
        dNa = d + 2*ha
        dan = dn + dNa - d
        den = 2*NOT/(abs(NOT))* math.sqrt( (math.sqrt((dan/2)**2 - (dbn/2)**2) - (math.pi*d*math.cos(helix_angle*math.pi/180)*math.cos(alpha_n*math.pi/180))/abs(NOT) * (eps_alphan-1))**2 + (dbn/2)**2) 
        alpha_en = math.acos(dbn/den)
        gamma_e = (0.5*math.pi+2*prof_shift_coeff_generated*math.tan(alpha_n*math.pi/180))/Zn + math.tan(alpha_n*math.pi/180) - alpha_n*math.pi/180 - math.tan(alpha_en) + alpha_en
        alpha_Fen = alpha_en - gamma_e
        hFe = (math.cos(gamma_e)-math.sin(gamma_e)*math.tan(alpha_Fen))*den/2 -Y
        if overlap_ratio == 0 and eps_alphan < 2:
            f_eps = 1
        if overlap_ratio == 0 and eps_alphan >= 2:
            f_eps = 0.7
        if 0<overlap_ratio<1 and eps_alphan < 2:
            f_eps = (1 - overlap_ratio + overlap_ratio/eps_alphan)**(1/2)
        if 0 < overlap_ratio < 1 and eps_alphan >=2:
            f_eps = ((1-overlap_ratio)/2 + overlap_ratio/eps_alphan)**(1/2)
        if overlap_ratio>=1:
            f_eps = eps_alphan**(-0.5)
        Y_F = (6*hFe/m_n * math.cos(alpha_Fen))/((sFn/m_n)**2 * math.cos(alpha_n*math.pi/180))*f_eps
    return Y_F, rho_F, sFn, hFe
   
def calculate_Y_S_ISO(sFn, rho_F, hFe): #according to ISO 6336:3 2019
   L = sFn/hFe
   qs = sFn/(2*rho_F)
   return abs((1.2 + 0.13*L)*qs**(1/(1.21+(2.3/L))))

def calculate_Y_B_ISO(Kind_of_gear, m_n, sR): #according to ISO 6336:3 2019
    ht = 2.25*m_n
    if Kind_of_gear =='int':
        if sR/m_n >= 3.5:
            Y_B = 1.0
        elif 1.75 < sR/m_n < 3.5:
            Y_B = 1.15* math.log(8.324*m_n/sR)
        else:
            print('Error, this is not a good configuration according to ISO 6336, sR/m_n should be higher.')
    else:
        if sR/ht >=1.2:
            Y_B = 1.0
        elif 0.5 < sR/ht < 1.2: 
            Y_B = 1.6* math.log(2.242*ht/sR)
        else:
           print('Error, this is not a good configuration according to ISO 6336, sR/ht should be higher.')
    return Y_B

def Y_delta_rel_T_ISO(qs, rho_, Y_S, speed): #according to ISO 6336:3 2019
    chi = 1/5 * (1+2*qs)
    chi__ = 1/5*(1+2*2.5)
    Y_delta_rel_T = (1+np.sqrt(rho_*chi))/(1+np.sqrt(rho_*chi__))
    if speed == 0: #different for static stress
        Y_delta_rel_T = 0.44*Y_S + 0.12
    return Y_delta_rel_T
  
def Y_R_rel_T_iso(number_of_cycles, material): #some kind of interpolation for the Y_R_rel_T factor since it depends on the speed in KISSsoft
    if material == 'Steel 40CrMoV13-9':
        points = [(0.5e6, 0.992), (1e6, 0.991), (3e6, 0.990), (5e6, 0.990)]
    else:
        points = [(0.2e6, 0.971), (0.5e6, 0.966), (1e6, 0.963), (3e6, 0.957), (5e6, 0.957)]

    # Sort points in case they are not in order
    points.sort()

    # Clamp if outside the known range
    if number_of_cycles <= points[0][0]:
        return points[0][1]
    elif number_of_cycles >= points[-1][0]:
        return points[-1][1]

    # Interpolate between the closest two points
    for i in range(1, len(points)):
        x0, y0 = points[i - 1]
        x1, y1 = points[i]
        if x0 <= number_of_cycles <= x1:
            # Linear interpolation
            return y0 + (y1 - y0) * (number_of_cycles - x0) / (x1 - x0)



def calculate_Y_NT_ISO(number_of_cycles): #interpolation accoding to ISO 6336:3 2019
    points = [
        (1e3, 2.5),
        (3e6, 1.0),
        (1e10, 0.85)
    ]

    # If below the lowest point, return highest value
    if number_of_cycles <= points[0][0]:
        return points[0][1]
    
    # If above the highest point, return lowest value
    if number_of_cycles >= points[-1][0]:
        return points[-1][1]

    # Log-log interpolation
    for i in range(len(points) - 1):
        x0, y0 = points[i]
        x1, y1 = points[i + 1]
        if x0 <= number_of_cycles <= x1:
            log_x = np.log10(number_of_cycles)
            log_x0 = np.log10(x0)
            log_x1 = np.log10(x1)
            log_y0 = np.log10(y0)
            log_y1 = np.log10(y1)
            
            log_y = log_y0 + (log_y1 - log_y0) * ((log_x - log_x0) / (log_x1 - log_x0))
            return 10 ** log_y

def calculate_ZH_ISO(helix_angle, working_pres_angle, alpha_n, type_of_gears): #according to ISO 6336:2 2019
    alpha_t = np.arctan(np.tan(alpha_n*np.pi/180)/np.cos(helix_angle*np.pi/180))
    beta_b = np.arctan(np.tan(helix_angle*np.pi/180)*np.cos(alpha_n*np.pi/180))
    if type_of_gears == 'ext':
        beta_b = np.arcsin(np.sin(helix_angle*np.pi/180)*np.cos(alpha_n*np.pi/180))
    return np.sqrt(2*np.cos(beta_b)*np.cos(working_pres_angle)/(np.cos(alpha_t)**2*np.sin(working_pres_angle)))

def calculate_Z_B_ISO(working_pressure_angle_t, helix_angle, z1, z2,m_t, x_1, x_2, Epsilon_alpha, type_of_gears): #Z1 and z2 positive for external, neg for internal, according to ISO 6336:2 2019
    if type_of_gears == 'int':
        ha1 = (1+x_1)*m_t*np.cos(helix_angle*np.pi/180)
        ha2 = (1+x_2)*m_t*np.cos(helix_angle*np.pi/180)
        dNa1 = m_t*z1 + 2*ha1 #tip diameters
        dNa2 = m_t*-z2 - 2*ha2
        db1 = m_t*z1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(helix_angle*np.pi/180))) #base diameter = d*cos(alpha_t)
        db2 = m_t*-z2*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(helix_angle*np.pi/180)))
    else:
        ha1 = (1+x_1)*m_t*np.cos(helix_angle*np.pi/180)
        ha2 = (1+x_2)*m_t*np.cos(helix_angle*np.pi/180)
        dNa1 = m_t*z1 + 2*ha1 #tip diameters
        dNa2 = m_t*z2 + 2*ha2
        db1 = m_t*z1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(helix_angle*np.pi/180))) #base diameter = d*cos(alpha_t)
        db2 = m_t*z2*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(helix_angle*np.pi/180)))
    M1 = math.tan(working_pressure_angle_t*np.pi/180)/(np.sqrt( (np.sqrt(dNa1**2/db1**2 - 1) - 2*np.pi/z1)*(np.sqrt(dNa2**2/db2**2 - 1) - (Epsilon_alpha-1)*2*np.pi/z2)))
    if 2.1 < Epsilon_alpha <2.5:
        Z_B = 1 
    elif Epsilon_alpha>1 and helix_angle == 0:
        if M1 <= 1:
            Z_B = 1
        else:
            Z_B = M1
    elif helix_angle != 0 and overlap_ratio >= 1 and Epsilon_alpha>1:
        Z_B = np.sqrt(1.2)
    elif helix_angle !=0 and Epsilon_alpha >1 and overlap_ratio <1:
        if M1 <= 1:
            Z_B = 1+overlap_ratio*(np.sqrt(1.2)-1)
        else:
            Z_B = M1 + overlap_ratio*(np.sqrt(1.2)-M1)
    return Z_B

def calculate_Z_D_ISO(working_pressure_angle_t, helix_angle, z1, z2,m_t, x, Epsilon_alpha): #Z1 and z2 positive for external, neg for internal, according to ISO 6336:2 2019
    ha = (1+x)*m_t*np.cos(helix_angle*np.pi/180)
    dNa1 = m_t*z1 + 2*ha
    dNa2 = m_t*z2 + 2*ha
    db1 = m_t*z1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(helix_angle*np.pi/180))) #base diameter = d*cos(pressure_Angle)
    db2 = m_t*z2*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(helix_angle*np.pi/180)))
    M2 = math.tan(working_pressure_angle_t)/(np.sqrt( (np.sqrt(dNa2**2/db2**2 - 1) - 2*np.pi/z2)*(np.sqrt(dNa1**2/db1**2 - 1) - (Epsilon_alpha-1)*2*np.pi/z1)))
    if Epsilon_alpha>1 and helix_angle == 0:
        if module2_ok <= 1:
            Z_D = 1 
        else:
            Z_D = M2
    if helix_angle != 0 and overlap_ratio >= 1 and Epsilon_alpha>1:
        Z_D = np.sqrt(1.2)
    if helix_angle !=0 and Epsilon_alpha >1 and overlap_ratio <1:
        if M2 <= 1:
            Z_D = 1+overlap_ratio*(np.sqrt(1.2)-1)
        else:
            Z_D = M2 + overlap_ratio*(np.sqrt(1.2)-M2)
    return Z_D

def calculate_ZE_ISO(Material_properties , material): #according to ISO 6336:2 2019
    E_mod = Material_properties[material][0]
    poisson_ratio = Material_properties[material][1]
    return np.sqrt(E_mod/(2*np.pi*(1-poisson_ratio**2)))

def calculate_Z_eps_ISO(helix_angle, Epsilon_alpha): #according to ISO 6336:2 2019
    if helix_angle == 0:
        Z_eps = np.sqrt((4-Epsilon_alpha)/3)
    else:
        if overlap_ratio <1:
            Z_eps = np.sqrt((4-Epsilon_alpha)/3*(1-overlap_ratio) + overlap_ratio/Epsilon_alpha)
        else:
            Z_eps = np.sqrt(1/Epsilon_alpha)
    return Z_eps

def calculate_Z_beta_ISO(helix_angle): #according to ISO 6336:2 2019
    return np.sqrt(1/np.cos(helix_angle*np.pi/180))

def calculate_ZL_ISO(material): #according to ISO 6336:2 2019
    global C_Z_L
    v40 = 250 #mm²/s is the kinematic viscosity of the  Kluberfluid C-F 1 Ultra
    if 850 <= Material_properties[material][7] <= 1200:
        C_Z_L = Material_properties[material][7] / 437.5 + 0.6357
    elif 850 > Material_properties[material][7]:
        C_Z_L = 0.83
    elif Material_properties[material][7] > 1200:
        C_Z_L = 0.91
    Z_L = C_Z_L + 4*(1-C_Z_L)/(1.2+134/v40)**2
    return Z_L

def calculate_Z_V_ISO(vw): #according to ISO 6336:2 2019
    C_Zv = C_Z_L + 0.02
    ZV = C_Zv + (2*(1-C_Zv)/np.sqrt(0.8 + 32/vw))
    return ZV

def calculate_Z_NT_ISO(number_of_cycles): #same function as Y_NT but for Z_NT, according to ISO 6336:2 2019
    # Reference points (log-log scale)
    points = [
        (1e5, 1.6),
        (5e7, 1.0),
        (1e10, 0.85)
    ]

    # If below the lowest point, return highest value
    if number_of_cycles <= points[0][0]:
        return points[0][1]
    
    # If above the highest point, return lowest value
    if number_of_cycles >= points[-1][0]:
        return points[-1][1]

    # Log-log interpolation
    for i in range(len(points) - 1):
        x0, y0 = points[i]
        x1, y1 = points[i + 1]
        if x0 <= number_of_cycles <= x1:
            log_x = np.log10(number_of_cycles)
            log_x0 = np.log10(x0)
            log_x1 = np.log10(x1)
            log_y0 = np.log10(y0)
            log_y1 = np.log10(y1)
            
            log_y = log_y0 + (log_y1 - log_y0) * ((log_x - log_x0) / (log_x1 - log_x0))
            return 10 ** log_y

def calculate_Z_R_ISO(Z_1,Z_2, material, m_t_pinion, m_t_wheel, alpha_t, alpha_wt): #according to ISO 6336:2-2019
    RZ = Material_properties[material][8]
    #db must be negative for internal gears, pos for externel gears
    rho_1 = 0.5*m_t_pinion*Z_1*np.cos(alpha_t)*np.tan(alpha_wt)
    rho_2 = 0.5*m_t_wheel*Z_2*np.cos(alpha_t)*np.tan(alpha_wt)
    rho_red = (rho_1*rho_2)/(rho_1+rho_2)
    RZ10 = RZ*np.cbrt(10/rho_red*10**(-3))
    
    if 850 <= Material_properties[material][7] <= 1200:
        CZR = 0.32 - 0.0002*Material_properties[material][7] 
    elif 850 > Material_properties[material][7]:
        CZR = 0.15
    else:
        CZR = 0.08
    Z_R = (3/(RZ10))**CZR
    return Z_R


def show_configuration_number():
    # Display the "Chosen configuration number" input below the "Helix angle output stage" inputs when the solution for the basic wolfrom part are given
    config_label.grid(row=8, column=2, padx=10, pady=10, sticky="e")
    config_entry.grid(row=8, column=3, padx=10, pady=10, sticky="w")

def show_configuration_number_pre():
    # Display the "Chosen configuration number" input below the "Helix angle output stage" inputs when the solutions for the pregearing are given
    config_label_pre.grid(row=8, column=4, padx=10, pady=10, sticky="e")
    config_entry_pre.grid(row=8, column=5, padx=10, pady=10, sticky="w")
    button3.grid_forget()
    button4.grid(row=10, column=4, columnspan=2, padx=10, pady=10, sticky="s")

def create_gui(): #this function creates the GUI for the Wolfrom gearbox calculation program and sets up the layout and input fields
    global entry_boxes, window, material_var, hole_var, hole_diameter_entry, hole_diameter_label, config_label, config_entry, gear_var, config_label_pre, config_entry_pre, button4, allowance_var, button3

    window = tk.Tk()
    window.title('Python program to calculate a Wolfrom based gearbox')
    window.geometry("975x550+10+150") 

    # Define inputs by category for each column
    general_inputs = ["Max Pitch diameter (mm)", "Output torque (Nm)", 
                      "Output speed (rpm)"]
    wolfrom_inputs = ["Axial length (mm)","Desired wolfrom gear ratio", "Percentage deviation (%)", "Number of planets", 
                      "Module in input stage (mm)", "Module in output stage (mm)", "Minimum efficiency (%)"]
    pre_gearing_inputs = ["Axial length (mm)","Desired gear ratio pregearing", "Percentage deviation (%)", "Number of planets", 
                          "Module in input stage (mm)", "Module in output stage (mm)", "Minimum efficiency (%)"]

    entry_boxes = []

    # Configure rows and columns for equal distribution
    for i in range(11):  # Total rows (max inputs + labels + button)
        window.grid_rowconfigure(i, weight=1)

    for i in range(3):  # Total columns
        window.grid_columnconfigure(i * 2, weight=1)
        window.grid_columnconfigure(i * 2 + 1, weight=1)

    # Column titles
    title_labels = ['General inputs', 'Basic wolfrom inputs', 'Pre-gearing inputs']
    for i, title in enumerate(title_labels):
        title_label = tk.Label(window, text=title, font=("Arial", 14, "bold"))
        title_label.grid(row=0, column=i * 2, columnspan=2, padx=10, pady=10, sticky="n")

    # General inputs
    for i, label_text in enumerate(general_inputs):
        label = tk.Label(window, text=label_text)
        label.grid(row=i + 1, column=0, padx=10, pady=10, sticky="e")
        
        entry = tk.Entry(window)
        entry.grid(row=i + 1, column=1, padx=10, pady=10, sticky="w")
        entry_boxes.append(entry)

    allowance_label = tk.Label(window, text="Tooth thickness allowances")
    allowance_label.grid(row=4, column=0, padx=10, pady=10, sticky="e")

    allowance_var = tk.StringVar(window)
    allowances = ["DIN 3967 cd25", "DIN 3967 cd26", "DIN 3967 f22"]
    allowance_var.set(allowances[0])

    allowance_dropdown = tk.OptionMenu(window, allowance_var, *allowances)
    allowance_dropdown.grid(row=4,column=1, padx=10, pady=10, sticky="w")

    material_label = tk.Label(window, text="Material")
    material_label.grid(row=5, column=0, padx=10, pady=10, sticky="e")

    material_var = tk.StringVar(window)
    materials = ["Steel 20MNCr-5", "Steel 18CrNiMo7-6", "Steel 16MnCrS5", "Steel 40CrMoV13-9"]
    material_var.set(materials[0])

    material_dropdown = tk.OptionMenu(window, material_var, *materials)
    material_dropdown.grid(row=5,column=1, padx=10, pady=10, sticky="w")

    gear_label = tk.Label(window, text="Type of gear")
    gear_label.grid(row=6, column=0, padx=10, pady=10, sticky="e")

    gear_var = tk.StringVar(window)  # Default selection
    gear_var.set('Spur')

    gear_dropdown =  tk.OptionMenu(window, gear_var, "Spur", "Helical, eps_beta = 1", "Helical, eps_beta = 2")
    gear_dropdown.grid(row = 6, column = 1, padx = 10, pady = 10, sticky = "w")

    hole_label = tk.Label(window, text="Cable hole in middle?")
    hole_label.grid(row=7, column=0, padx=10, pady=10, sticky="e")

    hole_var = tk.StringVar(window)
    hole_var.set("No")

    hole_dropdown = tk.OptionMenu(window, hole_var, "Yes", "No", command=show_hole_input)
    hole_dropdown.grid(row=7, column=1, padx=10, pady=10, sticky="w")

    # Hole diameter input (initially hidden)
    hole_diameter_label = tk.Label(window, text="Hole diameter (mm)")
    hole_diameter_entry = tk.Entry(window)

    # Wolfrom inputs
    for i, label_text in enumerate(wolfrom_inputs):
        label = tk.Label(window, text=label_text)
        label.grid(row=i + 1, column=2, padx=10, pady=10, sticky="e")
        
        entry = tk.Entry(window)
        entry.grid(row=i + 1, column=3, padx=10, pady=10, sticky="w")
        entry_boxes.append(entry)

    # Pre-gearing inputs
    for i, label_text in enumerate(pre_gearing_inputs):
        label = tk.Label(window, text=label_text)
        label.grid(row=i + 1, column=4, padx=10, pady=10, sticky="e")
        
        entry = tk.Entry(window)
        entry.grid(row=i + 1, column=5, padx=10, pady=10, sticky="w")
        entry_boxes.append(entry)

    # Configuration number label and entry (initially hidden)
    config_label = tk.Label(window, text="Basic Wolfrom configuration")
    config_entry = tk.Entry(window)
    entry_boxes.append(config_entry)

    config_label_pre = tk.Label(window, text="Pregearing configuration")
    config_entry_pre = tk.Entry(window)
    entry_boxes.append(config_entry_pre)

    # Buttons in the second and third columns
    button1 = tk.Button(window, text="Preliminary study", command=calculate_preliminary)
    button1.grid(row=10, column=0, columnspan=2, padx=10, pady=10, sticky="s")

    button2 = tk.Button(window, text="Calculate Wolfrom Gearbox", command=calculate_basic_wolfrom)
    button2.grid(row=10, column=2, columnspan=2, padx=10, pady=10, sticky="s")

    button3 = tk.Button(window, text="Calculate Pre-Gearing", command=calculate_pre_gearing)
    button3.grid(row=10, column=4, columnspan=2, padx=10, pady=10, sticky="s")

    button4 = tk.Button(window, text="Summary", command= summary)

    window.mainloop()

def show_hole_input(selection): # Function to show or hide the hole diameter input based on the selection
    global hole_diameter_entry, hole_diameter_label
    if selection == "Yes":
        # Create label and entry for hole diameter if "Yes" is selected
        hole_diameter_label = tk.Label(window, text="Hole diameter (mm)")
        hole_diameter_label.grid(row=8, column=0, padx=10, pady=10, sticky="e")
        
        hole_diameter_entry = tk.Entry(window)
        hole_diameter_entry.grid(row=8, column=1, padx=10, pady=10, sticky="w")
        entry_boxes.append(hole_diameter_entry)
    else:
        # Remove hole diameter entry if "No" is selected (if exists)
        if 'hole_diameter_entry' in globals():
            hole_diameter_entry.grid_forget()
            hole_diameter_label.grid_forget()

def calculate_width(Ft1,Ft2, axial_length, type_of_gear, mt1, mt2): #a function to calculate the width of each stage, done differently for spur and helical gears
    if type_of_gear == 'Spur':
        Ft1 = abs(Ft1)
        Ft2 = abs(Ft2)
        b2 = axial_length/(Ft1/Ft2 + 1)
        b1 = (Ft1*b2)/(Ft2)
    else:
        Ft1 = abs(Ft1)
        Ft2 = abs(Ft2)
        b2 = (Ft2*mt2*axial_length/(Ft1*mt1+ Ft2*mt2))
        b1 = axial_length-b2
    return b1, b2 #in mm

def calculate_F_a(F_t, helix_angle): #calculate the axial force from the tangential force, according to the ISO 6336 norm
    return F_t*np.tan(helix_angle*np.pi/180)

def calculate_F_r(F_t, pressure_angle, helix_angle): #calculate the radial force from the tangential force, according to the ISO 6336 norm
    return F_t*np.tan(pressure_angle*np.pi/180)/np.cos(helix_angle*np.pi/180)

def calculate_F_n(F_t, pressure_angle): #calculate the normal force from the tangential force, according to the ISO 6336 norm
    return F_t/np.cos(pressure_angle*np.pi/180)

#profile shift function which uses the decrease/increase in center distance, based on the neighboring condition
def profile_shift_basic_wolfrom(material, Possible_solutions_wolfrom, m_t_1, m_t_2, NOP, allowances):
    possibilities = [] #some initial lists, used to store solutions
    possibilities_2 = []
    bending_stresses_after_PS = []
    Z_E = calculate_ZE_ISO(Material_properties, material)
    for i in range(len(Possible_solutions_wolfrom)):
        exit_inner_loops = False
        exit_inner_loops_second_time = False
        #calculation of the base diameters for each gear
        db2_2 = m_t_2*10**(-3)*Possible_solutions_wolfrom[i][0]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) #for the internal one => must be alpha_t!
        db1_2 = m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) #for the external one => must be alpha_t!
        db2_1= m_t_1*10**(-3)*Possible_solutions_wolfrom[i][1]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))) #=> must be alpha_t!
        db1_1 = m_t_1*10**(-3)*Possible_solutions_wolfrom[i][2]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))) #=> must be alpha_t!
        tau_1 = (2*np.pi)/(Possible_solutions_wolfrom[i][2])
        tau_2 = (2*np.pi)/(Possible_solutions_wolfrom[i][3])
        y_max_1 = (pitch_diam*10**(-3)/2)/(m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) - (Possible_solutions_wolfrom[i][2]*m_t_1*10**(-3))/(2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) - (Possible_solutions_wolfrom[i][1] - Possible_solutions_wolfrom[i][2])/(2*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))  #maximum positive y-factor
        cos_alphawt_1 = np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)))/((2*y_max_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)/(Possible_solutions_wolfrom[i][1]-Possible_solutions_wolfrom[i][2])) + 1)
        alpha_wt_1 = np.arccos(cos_alphawt_1)
        inv_alpha_wt_min_1 = np.tan(alpha_wt_1) - alpha_wt_1
        max_sum_of_prof_shift_1 = (inv_alpha_wt_min_1-np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))) + np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)))/(2*np.tan(pressure_angle*np.pi/180))*(Possible_solutions_wolfrom[i][1] - Possible_solutions_wolfrom[i][2])
        if max_sum_of_prof_shift_1 > 1: #if higher than 2, the program takes too long to go over all the possibilities, thus it is bounded.
            max_sum_of_prof_shift_1 = 1
        if abs(round(max_sum_of_prof_shift_1,10)) == 0: #when no room for center distance to increase, the only possibility is a negative PS for the ring and a positive on the planet, which will not increase the safety factor for both the gears, thus the initial configuration with no profile shift is shown.
            exit_inner_loops = True
            bending_stresses_after_PS.append([Possible_solutions_wolfrom[i][37], Possible_solutions_wolfrom[i][38], Possible_solutions_wolfrom[i][26], Possible_solutions_wolfrom[i][25], 0, 0 , Possible_solutions_wolfrom[i][13],Possible_solutions_wolfrom[i][29], Possible_solutions_wolfrom[i][30], Possible_solutions_wolfrom[i][45], Possible_solutions_wolfrom[i][46]])
        x_P1 = 0
        sum = 0
        step = 0.005
        while sum <= abs(max_sum_of_prof_shift_1):
            if exit_inner_loops == False:
                while x_P1 >= -max_sum_of_prof_shift_1:
                    if x_P1 <1-(Possible_solutions_wolfrom[i][2]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)): #checking undercut
                        x_P1 -= step #if met the current division of PS is not good, so the next iteration will be done
                        if x_P1 < -max_sum_of_prof_shift_1: #if lower than this value, the loop is at its end
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1: #if outer loop also at its end
                                if len(possibilities)==0: #then there is no solution found, a list of zeros is returned, the code might need a rerun with a bigger difference allowed in the safety factors down below
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0])) #the solution with the highest value of the safety factor in the ring gear is chosen
                                    bending_stresses_after_PS.append(Best_config) #stored in a list
                                    possibilities.clear()
                                    exit_inner_loops = True #time to go to the output stage!
                                    exit_inner_loops_second_time = False
                                break
                        continue
                    x_R1 = -sum - x_P1
                    XE_1_ext = x_P1 + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) #computing generating profile shifts
                    XE_1_int = x_R1 + Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    inv_alphawt = 2*np.tan(pressure_angle*np.pi/180)*((-x_P1-x_R1)/(Possible_solutions_wolfrom[i][1]-Possible_solutions_wolfrom[i][2])) + np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    if inv_alphawt < 0: #because then the inverse involute is not converging and it also does not make sense to have a negative value
                        x_P1 -= step
                        if x_P1 < -max_sum_of_prof_shift_1:
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1:
                                if len(possibilities)==0:
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0]))
                                    bending_stresses_after_PS.append(Best_config)
                                    possibilities.clear()
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                break
                        continue
                    alphawt = inverse_involute(inv_alphawt)
                    y = (Possible_solutions_wolfrom[i][1] - Possible_solutions_wolfrom[i][2])/(2*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))*((np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)))/np.cos(alphawt))-1)
                    da2_1 = m_t_1*10**(-3)*Possible_solutions_wolfrom[i][1] - 2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*(1+x_R1) #tip diameters
                    da1_1 = m_t_1*10**(-3)*Possible_solutions_wolfrom[i][2] + 2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*(1+x_P1)
                    a_1 = ((Possible_solutions_wolfrom[i][1]-Possible_solutions_wolfrom[i][2])/(2*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) + y)*m_t_1*10**-3*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)
                    alpha_ta1 = np.arccos(db1_1/da1_1) #planet tangential pressure angle at the tip
                    alpha_ta2 = np.arccos(db2_1/da2_1) #ring tangential pressure angle at the tip
                    beta_a_planeet = np.arctan((np.tan(Possible_solutions_wolfrom[i][9]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))))/np.cos(alpha_ta1))
                    beta_a_ring = np.arctan((np.tan(Possible_solutions_wolfrom[i][9]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))))/np.cos(alpha_ta2))
                    if da1_1*(np.pi/(2*Possible_solutions_wolfrom[i][2]) + 2*x_P1*np.tan(pressure_angle*np.pi/180)/Possible_solutions_wolfrom[i][2] + np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) - np.tan(alpha_ta1)+ alpha_ta1)*np.cos(beta_a_planeet) < 0.4*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180): #tooth thickness check for planet
                        x_P1 -= step
                        if x_P1 < -max_sum_of_prof_shift_1:
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1:
                                if len(possibilities)==0:
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0]))
                                    bending_stresses_after_PS.append(Best_config)
                                    possibilities.clear()
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                break
                        continue
                    if da2_1*(np.pi/(2*Possible_solutions_wolfrom[i][1]) + 2*x_R1*np.tan(pressure_angle*np.pi/180)/Possible_solutions_wolfrom[i][1] - (np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) - np.tan(alpha_ta2)+ alpha_ta2))*np.cos(beta_a_ring) < 0.4*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180): #tooth thickness check for ring
                        x_P1 -= step
                        if x_P1 < -max_sum_of_prof_shift_1:
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1:
                                if len(possibilities)==0:
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0]))
                                    bending_stresses_after_PS.append(Best_config)
                                    possibilities.clear()
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                break
                        continue
                    #calculation of the transverse contact ratio from the values of the profile shift coefficients, according to the ISO 6336 norm
                    d_root_ring =  m_t_1*10**(-3)*Possible_solutions_wolfrom[i][1] + 2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*(1.25 - x_R1) #calculation of root diameters for ring and planet
                    d_root_planet =  m_t_1*10**(-3)*Possible_solutions_wolfrom[i][2] - 2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*(1.25 - x_P1)
                    dw_p1 = db1_1/np.cos(alphawt) #calculation of working pitch circles for planet and ring
                    dw_ring1 = db2_1/np.cos(alphawt)
                    aw = 1/2*(dw_ring1 - dw_p1)
                    dNf1 = np.sqrt((2*aw*np.sin(alphawt) - np.sqrt(da2_1**2 - db2_1**2))**2 + db1_1**2) #active root diameter for planet (1) and ring (2)
                    dNf2 = np.sqrt((2*aw*np.sin(alphawt) + np.sqrt(da1_1**2 - db1_1**2))**2 + db2_1**2)
                    ksi_Nfw1 = np.tan(alphawt) - np.tan(np.arccos(db1_1/dNf1))
                    ksi_Nfw2 = np.tan(alphawt) - np.tan(np.arccos(db2_1/dNf2))
                    ksiNaw1 = ksi_Nfw2*(-Possible_solutions_wolfrom[i][1]/Possible_solutions_wolfrom[i][2])
                    Eps_alpha = (ksi_Nfw1 + ksiNaw1)/tau_1 
                    if a_1 + da1_1/2 > d_root_ring/2 - 0.2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180) or da2_1/2 < d_root_planet/2 + a_1 +0.2*m_t_1*10**(-3)*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180) or db2_1 >= da2_1:# interference check
                        x_P1 -= step
                        if x_P1 < -max_sum_of_prof_shift_1:
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1:
                                if len(possibilities)==0:
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0]))
                                    bending_stresses_after_PS.append(Best_config)
                                    possibilities.clear()
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                break
                        continue
                    if Eps_alpha < 1 or np.isnan(Eps_alpha): #norm doesn't say anything about transverse contact ratio's smaller than 1, so these solutions must be avoided
                        x_P1 -= step
                        if x_P1 < -max_sum_of_prof_shift_1:
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1:
                                if len(possibilities)==0:
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0]))
                                    bending_stresses_after_PS.append(Best_config)
                                    possibilities.clear()
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                break
                        continue
                    #calculation of all the bending stresses
                    Y_F_R1,rho_F_R1, sFn_R1, hFe_R1  = calculate_Y_F_ISO('int',Possible_solutions_wolfrom[i][1], Possible_solutions_wolfrom[i][9], XE_1_int, m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_1*10**-3*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180), Eps_alpha, x_R1)
                    Y_F_p1,rho_F_p1, sFn_p1, hFe_p1  = calculate_Y_F_ISO('ext',Possible_solutions_wolfrom[i][2], Possible_solutions_wolfrom[i][9], XE_1_ext, m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_1*10**-3*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180), Eps_alpha, x_P1)
                    sigma_R1 = Possible_solutions_wolfrom[i][21]/NOP/(Possible_solutions_wolfrom[i][7]*10**-3*m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*10**-3)*Y_F_R1*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_wolfrom[i][9])*calculate_Y_B_ISO('int', m_t_1*10**-3*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180), 3.5*m_t_1*10**-3*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R1, rho_F_R1, hFe_R1)*K
                    sigma_p1 = Possible_solutions_wolfrom[i][22]/(Possible_solutions_wolfrom[i][7]*10**-3*m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*10**-3)*Y_F_p1*calculate_Y_S_ISO(sFn_p1, rho_F_p1, hFe_p1)*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_wolfrom[i][9])*calculate_Y_B_ISO('ext', m_t_1*10**-3*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180), 3.5*m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)*10**-3)*Y_DT*K
                    sigma_FP_R1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP)*Y_delta_rel_T_ISO(sFn_R1/(2*rho_F_R1), Material_properties[material][5], calculate_Y_S_ISO(sFn_R1, rho_F_R1, hFe_R1),1)*Y_R_rel_T_iso(3e6*NOP, material)*Y_X/SF_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    sigma_FP_P1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p1/(2*rho_F_p1), Material_properties[material][5], calculate_Y_S_ISO(sFn_p1, rho_F_p1, hFe_p1),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    sigma_FG_R1 = sigma_FP_R1*SF_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    sigma_FG_p1 = sigma_FP_P1*SF_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    #calculation of all the Hertzian stresses
                    sigma_H_R1 = calculate_ZH_ISO(Possible_solutions_wolfrom[i][9], alphawt, pressure_angle, 'int')*Z_E*calculate_Z_eps_ISO(Possible_solutions_wolfrom[i][9], Eps_alpha)*calculate_Z_beta_ISO(Possible_solutions_wolfrom[i][9])*np.sqrt((Possible_solutions_wolfrom[i][21]/NOP)/(m_t_1*10**(-3)*Possible_solutions_wolfrom[i][2]*Possible_solutions_wolfrom[i][7]*10**(-3))*(-Possible_solutions_wolfrom[i][1]/Possible_solutions_wolfrom[i][2]+1)/(-Possible_solutions_wolfrom[i][1]/Possible_solutions_wolfrom[i][2]))*np.sqrt(K) #Z_D is = 1 for internal gears
                    sigma_H_P1 = sigma_H_R1*calculate_Z_B_ISO(alphawt*180/np.pi, Possible_solutions_wolfrom[i][9], Possible_solutions_wolfrom[i][2], -Possible_solutions_wolfrom[i][1],m_t_1*10**(-3), x_P1, x_R1, Eps_alpha, 'int')
                    sigma_HP_R1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP)/SH_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_wolfrom[i][2],-Possible_solutions_wolfrom[i][1], material, m_t_1*10**(-3), m_t_1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)), alphawt)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[i][6])
                    sigma_HP_P1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_wolfrom[i][2],-Possible_solutions_wolfrom[i][1], material, m_t_1*10**(-3), m_t_1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)), alphawt)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[i][6])
                    sigma_HG_R1 = sigma_HP_R1*SH_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    sigma_HG_p1 = sigma_HP_P1*SH_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180))
                    if sigma_FG_R1/(sigma_R1*10**(-6)) > SF_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)) and  sigma_FG_p1/(sigma_p1*10**(-6)) > SF_min(m_t_1*np.cos(Possible_solutions_wolfrom[i][9]*np.pi/180)): #check if the safety factors are higher than the minimum ones
                        if abs(sigma_FG_R1/(sigma_R1*10**(-6)) -  sigma_FG_p1/(sigma_p1*10**(-6)))<= 0.02: #'equalized' safety factors check
                            possibilities.append([sigma_FG_R1/(sigma_R1*10**(-6)), sigma_FG_p1/(sigma_p1*10**(-6)), sigma_p1*10**(-6), sigma_R1*10**(-6), x_P1, x_R1 , a_1, sigma_H_R1*10**(-6), sigma_H_P1*10**(-6), sigma_HG_R1/(sigma_H_R1*10**(-6)), sigma_HG_p1/(sigma_H_P1*10**(-6))]) #solution is stored
                            x_P1 -= step
                            if x_P1 < -max_sum_of_prof_shift_1:
                                sum += 0.01
                                x_P1 = 0
                                if sum > max_sum_of_prof_shift_1:
                                    if len(possibilities)==0:
                                        bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                        exit_inner_loops = True
                                        exit_inner_loops_second_time = False
                                    else:
                                        Best_config = max(possibilities, key=lambda x:(x[0]))
                                        bending_stresses_after_PS.append(Best_config)
                                        possibilities.clear()
                                        exit_inner_loops = True
                                        exit_inner_loops_second_time = False
                                    break
                            continue
                        else: #if not close enough to each other
                            x_P1 -= step
                            if x_P1 < -max_sum_of_prof_shift_1:
                                sum += 0.01
                                x_P1 = 0
                                if sum > max_sum_of_prof_shift_1:
                                    if len(possibilities)==0:
                                        bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                        exit_inner_loops = True
                                        exit_inner_loops_second_time = False
                                    else:
                                        Best_config = max(possibilities, key=lambda x:(x[0]))
                                        bending_stresses_after_PS.append(Best_config)
                                        possibilities.clear()
                                        exit_inner_loops = True
                                        exit_inner_loops_second_time = False
                                    break
                            continue
                    else: #if not strong enough, go to next iteration
                        x_P1 -= step
                        if x_P1 < -max_sum_of_prof_shift_1:
                            sum += 0.01
                            x_P1 = 0
                            if sum > max_sum_of_prof_shift_1:
                                if len(possibilities)==0:
                                    bending_stresses_after_PS.append([0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0])
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                else:
                                    Best_config = max(possibilities, key=lambda x:(x[0]))
                                    bending_stresses_after_PS.append(Best_config)
                                    possibilities.clear()
                                    exit_inner_loops = True
                                    exit_inner_loops_second_time = False
                                break
                        continue
            if exit_inner_loops:
                if bending_stresses_after_PS[i][4] == 0 and bending_stresses_after_PS[i][5] == 0: #if no profile shift in input stage, the reference center distance is the original one
                    a = Possible_solutions_wolfrom[i][13]
                else:
                    a = Best_config[6] #otherwise it is the new center distance
                y_2 = abs(a-Possible_solutions_wolfrom[i][14])/(m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                if a < Possible_solutions_wolfrom[i][14]:
                    y_2 = -y_2
                a_2 = ((Possible_solutions_wolfrom[i][0] - Possible_solutions_wolfrom[i][3])/(2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) + y_2)*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)
                cos_alphawt_2 = np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)))/((2*y_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)/(Possible_solutions_wolfrom[i][0]-Possible_solutions_wolfrom[i][3])) + 1)
                alpha_wt_2 = np.arccos(cos_alphawt_2)
                inv_alpha_wt_min_2 = np.tan(alpha_wt_2) - alpha_wt_2
                sum_of_prof_shift_2 = (inv_alpha_wt_min_2-np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) + np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)))/(2*np.tan(pressure_angle*np.pi/180))*(Possible_solutions_wolfrom[i][0] - Possible_solutions_wolfrom[i][3]) 
                underlimit = sum_of_prof_shift_2
                if abs(round(sum_of_prof_shift_2,10)) == 0: #if no adjustment is needed, PS coefficients are zero and the
                    bending_stresses_after_PS[i].extend([Possible_solutions_wolfrom[i][36], Possible_solutions_wolfrom[i][39], Possible_solutions_wolfrom[i][27], Possible_solutions_wolfrom[i][24], 0, 0 , Possible_solutions_wolfrom[i][14], Possible_solutions_wolfrom[i][28], Possible_solutions_wolfrom[i][31], Possible_solutions_wolfrom[i][32], Possible_solutions_wolfrom[i][35],Possible_solutions_wolfrom[i][44],Possible_solutions_wolfrom[i][47]])
                    exit_inner_loops_second_time = True
                    break
                x_P2 = -sum_of_prof_shift_2 #same logic as in first part
                step = 0.005
                if sum_of_prof_shift_2 < 0 and exit_inner_loops_second_time == False:
                    while x_P2 >= sum_of_prof_shift_2:
                        x_R2 = -sum_of_prof_shift_2 - x_P2
                        XE_2_ext = x_P2 + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        XE_2_int = x_R2 + Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        da2_2 = m_t_2*10**(-3)*Possible_solutions_wolfrom[i][0] - 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1+x_R2)
                        da1_2 = m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3] + 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1+x_P2)
                        d_root_ring =  m_t_2*10**(-3)*Possible_solutions_wolfrom[i][0] + 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1.25 - x_R2)
                        d_root_planet =  m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3] - 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1.25 - x_P2)
                        dw_p2 = db1_2/np.cos(alpha_wt_2)
                        dw_ring2 = db2_2/np.cos(alpha_wt_2)
                        aw = 1/2*(dw_ring2 - dw_p2)
                        alpha_ta1 = np.arccos(db1_2/da1_2) #planeet
                        alpha_ta2 = np.arccos(db2_2/da2_2) #ring
                        beta_a_planeet = np.arctan((np.tan(Possible_solutions_wolfrom[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))))/np.cos(alpha_ta1))
                        beta_a_ring = np.arctan((np.tan(Possible_solutions_wolfrom[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))))/np.cos(alpha_ta2))
                        dNf1_2 = np.sqrt((2*aw*np.sin(alpha_wt_2) - np.sqrt(da2_2**2 - db2_2**2))**2 + db1_2**2) #planeet
                        dNf2_2 = np.sqrt((2*aw*np.sin(alpha_wt_2) + np.sqrt(da1_2**2 - db1_2**2))**2 + db2_2**2) #ring
                        ksi_Nfw1_2 = np.tan(alpha_wt_2) - np.tan(np.arccos(db1_2/dNf1_2)) #planeet
                        ksi_Nfw2_2 = np.tan(alpha_wt_2) - np.tan(np.arccos(db2_2/dNf2_2)) #ring
                        ksiNaw1_2 = ksi_Nfw2_2*(-Possible_solutions_wolfrom[i][0]/Possible_solutions_wolfrom[i][3])
                        Eps_alpha_2 = (ksi_Nfw1_2 + ksiNaw1_2)/tau_2
                        if Eps_alpha_2 <1 or np.isnan(Eps_alpha_2):
                            x_P2 -= step
                            if x_P2 < sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        if a_2 + da1_2/2 > d_root_ring/2 - 0.2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180) or da2_2/2 < d_root_planet/2 + a_2 +0.2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180) or db2_2 >= da2_2:
                            x_P2 -= step
                            if x_P2 < sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        if da1_2*(np.pi/(2*Possible_solutions_wolfrom[i][3]) + 2*x_P2*np.tan(pressure_angle*np.pi/180)/Possible_solutions_wolfrom[i][3] + np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) - np.tan(alpha_ta1)+ alpha_ta1)*np.cos(beta_a_planeet) < 0.4*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180):
                            x_P2 -= step
                            if x_P2 < sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        if da2_2*(np.pi/(2*Possible_solutions_wolfrom[i][0]) + 2*x_R2*np.tan(pressure_angle*np.pi/180)/Possible_solutions_wolfrom[i][0] - (np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) - np.tan(alpha_ta2)+ alpha_ta2))*np.cos(beta_a_ring) < 0.4*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180):
                            x_P2 -= step
                            if x_P2 < sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        Y_F_R2,rho_F_R2, sFn_R2, hFe_R2  = calculate_Y_F_ISO('int',Possible_solutions_wolfrom[i][0], Possible_solutions_wolfrom[i][10], XE_2_int, m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), Eps_alpha_2, x_R2)
                        Y_F_p2,rho_F_p2, sFn_p2, hFe_p2  = calculate_Y_F_ISO('ext',Possible_solutions_wolfrom[i][3], Possible_solutions_wolfrom[i][10], XE_2_ext, m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), Eps_alpha_2, x_P2)
                        sigma_R2 = Possible_solutions_wolfrom[i][20]/NOP/(Possible_solutions_wolfrom[i][8]*10**-3*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3)*Y_F_R2*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_wolfrom[i][10])*calculate_Y_B_ISO('int', m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), 3.5*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R2, rho_F_R2, hFe_R2)*K
                        sigma_p2 = Possible_solutions_wolfrom[i][23]/(Possible_solutions_wolfrom[i][8]*10**-3*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3)*Y_F_p2*calculate_Y_S_ISO(sFn_p2, rho_F_p2, hFe_p2)*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_wolfrom[i][10])*calculate_Y_B_ISO('ext', m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), 3.5*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3)*Y_DT*K
                        sigma_FP_R2 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP)*Y_delta_rel_T_ISO(sFn_R2/(2*rho_F_R2), Material_properties[material][5], calculate_Y_S_ISO(sFn_R2, rho_F_R2, hFe_R2),1)*Y_R_rel_T_iso(3e6*NOP, material)*Y_X/SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_FP_P2 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p2/(2*rho_F_p2), Material_properties[material][5], calculate_Y_S_ISO(sFn_p2, rho_F_p2, hFe_p2),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_FG_R2 = sigma_FP_R2*SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_FG_p2 = sigma_FP_P2*SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_H_R2 = calculate_ZH_ISO(Possible_solutions_wolfrom[i][10], alpha_wt_2, pressure_angle, 'int')*Z_E*calculate_Z_eps_ISO(Possible_solutions_wolfrom[i][10], Eps_alpha_2)*calculate_Z_beta_ISO(Possible_solutions_wolfrom[i][10])*np.sqrt((Possible_solutions_wolfrom[i][20]/NOP)/(m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3]*Possible_solutions_wolfrom[i][8]*10**(-3))*(-Possible_solutions_wolfrom[i][0]/Possible_solutions_wolfrom[i][3]+1)/(-Possible_solutions_wolfrom[i][0]/Possible_solutions_wolfrom[i][3]))*np.sqrt(K) #Z_D is = 1 for internal gears
                        sigma_H_P2 = sigma_H_R2*calculate_Z_B_ISO(alpha_wt_2*180/np.pi, Possible_solutions_wolfrom[i][10], Possible_solutions_wolfrom[i][3], -Possible_solutions_wolfrom[i][0],m_t_2*10**(-3), x_P2, x_R2, Eps_alpha_2, 'int')
                        sigma_HP_R2 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP)/SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_wolfrom[i][3],-Possible_solutions_wolfrom[i][0], material, m_t_2*10**(-3), m_t_2*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)), alpha_wt_2)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[i][5])
                        sigma_HP_P2 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_wolfrom[i][3],-Possible_solutions_wolfrom[i][0], material, m_t_2*10**(-3), m_t_2*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)), alpha_wt_2)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[i][5])
                        sigma_HG_R2 = sigma_HP_R2*SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_HG_p2 = sigma_HP_P2*SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        if sigma_FG_R2/(sigma_R2*10**(-6)) > SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) and  sigma_FG_p2/(sigma_p2*10**(-6)) > SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) and  sigma_FG_R2/(sigma_R2*10**(-6))>Possible_solutions_wolfrom[i][36] and sigma_FG_p2/(sigma_p2*10**(-6)) > Possible_solutions_wolfrom[i][39] :
                            possibilities_2.append([(sigma_FG_R2/(sigma_R2*10**(-6))), sigma_FG_p2/(sigma_p2*10**(-6)), sigma_p2*10**(-6), sigma_R2*10**(-6), x_P2, x_R2, a_2, sigma_H_R2*10**(-6), sigma_H_P2*10**(-6), sigma_HP_R2, sigma_HP_P2, sigma_HG_R2/(sigma_H_R2*10**(-6)), sigma_HG_p2/(sigma_H_P2*10**(-6))])
                            x_P2 -= step
                            if x_P2 < sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        else:
                            x_P2 -= step
                            if x_P2 < sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                if sum_of_prof_shift_2 >= 0 and exit_inner_loops_second_time == False:
                    x_P2 = sum_of_prof_shift_2
                    step = -0.005
                    underlimit = -sum_of_prof_shift_2
                    if -sum_of_prof_shift_2 < 1-(Possible_solutions_wolfrom[i][3]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)):
                        underlimit =  1-(Possible_solutions_wolfrom[i][3]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                    while x_P2 >= underlimit:
                        x_R2 = -sum_of_prof_shift_2 - x_P2
                        XE_2_ext = x_P2 + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        XE_2_int = x_R2 + Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        da2_2 = m_t_2*10**(-3)*Possible_solutions_wolfrom[i][0] - 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1+x_R2)
                        da1_2 = m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3] + 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1+x_P2)
                        d_root_ring =  m_t_2*10**(-3)*Possible_solutions_wolfrom[i][0] + 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1.25 - x_R2)
                        d_root_planet =  m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3] - 2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*(1.25 - x_P2)
                        dw_p2 = db1_2/np.cos(alpha_wt_2)
                        dw_ring2 = db2_2/np.cos(alpha_wt_2)
                        aw = 1/2*(dw_ring2 - dw_p2)
                        alpha_ta1 = np.arccos(db1_2/da1_2) #planeet
                        alpha_ta2 = np.arccos(db2_2/da2_2) #ring
                        beta_a_planeet = np.arctan((np.tan(Possible_solutions_wolfrom[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))))/np.cos(alpha_ta1))
                        beta_a_ring = np.arctan((np.tan(Possible_solutions_wolfrom[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))))/np.cos(alpha_ta2))
                        dNf1_2 = np.sqrt((2*aw*np.sin(alpha_wt_2) - np.sqrt(da2_2**2 - db2_2**2))**2 + db1_2**2) #planeet
                        dNf2_2 = np.sqrt((2*aw*np.sin(alpha_wt_2) + np.sqrt(da1_2**2 - db1_2**2))**2 + db2_2**2) #ring
                        ksi_Nfw1_2 = np.tan(alpha_wt_2) - np.tan(np.arccos(db1_2/dNf1_2)) #planeet
                        ksi_Nfw2_2 = np.tan(alpha_wt_2) - np.tan(np.arccos(db2_2/dNf2_2)) #ring
                        ksiNaw1_2 = ksi_Nfw2_2*(-Possible_solutions_wolfrom[i][0]/Possible_solutions_wolfrom[i][3])
                        Eps_alpha_2 = (ksi_Nfw1_2 + ksiNaw1_2)/tau_2
                        if Eps_alpha_2 <1 or np.isnan(Eps_alpha_2):
                            x_P2 += step
                            if x_P2 < underlimit:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        if a_2 + da1_2/2 > d_root_ring/2 - 0.2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180) or da2_2/2 < d_root_planet/2 + a_2 +0.2*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180) or db2_2 >= da2_2:
                            x_P2 += step
                            if x_P2 < -sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        if da1_2*(np.pi/(2*Possible_solutions_wolfrom[i][3]) + 2*x_P2*np.tan(pressure_angle*np.pi/180)/Possible_solutions_wolfrom[i][3] + np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) - np.tan(alpha_ta1)+ alpha_ta1)*np.cos(beta_a_planeet) < 0.4*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180):
                            x_P2 += step
                            if x_P2 < -sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        if da2_2*(np.pi/(2*Possible_solutions_wolfrom[i][0]) + 2*x_R2*np.tan(pressure_angle*np.pi/180)/Possible_solutions_wolfrom[i][0] - (np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) - np.tan(alpha_ta2)+ alpha_ta2))*np.cos(beta_a_ring) < 0.4*m_t_2*10**(-3)*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180):
                            x_P2 += step
                            if x_P2 < -sum_of_prof_shift_2:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                            continue
                        Y_F_R2,rho_F_R2, sFn_R2, hFe_R2  = calculate_Y_F_ISO('int',Possible_solutions_wolfrom[i][0], Possible_solutions_wolfrom[i][10], XE_2_int, m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), Eps_alpha_2, x_R2)
                        Y_F_p2,rho_F_p2, sFn_p2, hFe_p2  = calculate_Y_F_ISO('ext',Possible_solutions_wolfrom[i][3], Possible_solutions_wolfrom[i][10], XE_2_ext, m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), Eps_alpha_2, x_P2)
                        sigma_R2 = Possible_solutions_wolfrom[i][20]/NOP/(Possible_solutions_wolfrom[i][8]*10**-3*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3)*Y_F_R2*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_wolfrom[i][10])*calculate_Y_B_ISO('int', m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), 3.5*m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R2, rho_F_R2, hFe_R2)*K
                        sigma_p2 = Possible_solutions_wolfrom[i][23]/(Possible_solutions_wolfrom[i][8]*10**-3*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3)*Y_F_p2*calculate_Y_S_ISO(sFn_p2, rho_F_p2, hFe_p2)*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_wolfrom[i][10])*calculate_Y_B_ISO('ext', m_t_2*10**-3*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180), 3.5*m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)*10**-3)*Y_DT*K
                        sigma_FP_R2 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP)*Y_delta_rel_T_ISO(sFn_R2/(2*rho_F_R2), Material_properties[material][5], calculate_Y_S_ISO(sFn_R2, rho_F_R2, hFe_R2),1)*Y_R_rel_T_iso(3e6**NOP, material)*Y_X/SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_FP_P2 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p2/(2*rho_F_p2), Material_properties[material][5], calculate_Y_S_ISO(sFn_p2, rho_F_p2, hFe_p2),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_FG_R2 = sigma_FP_R2*SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_FG_p2 = sigma_FP_P2*SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))

                        sigma_H_R2 = calculate_ZH_ISO(Possible_solutions_wolfrom[i][10], alpha_wt_2, pressure_angle, 'int')*Z_E*calculate_Z_eps_ISO(Possible_solutions_wolfrom[i][10], Eps_alpha_2)*calculate_Z_beta_ISO(Possible_solutions_wolfrom[i][10])*np.sqrt((Possible_solutions_wolfrom[i][20]/NOP)/(m_t_2*10**(-3)*Possible_solutions_wolfrom[i][3]*Possible_solutions_wolfrom[i][8]*10**(-3))*(-Possible_solutions_wolfrom[i][0]/Possible_solutions_wolfrom[i][3]+1)/(-Possible_solutions_wolfrom[i][0]/Possible_solutions_wolfrom[i][3]))*np.sqrt(K) #Z_D is = 1 for internal gears
                        sigma_H_P2 = sigma_H_R2*calculate_Z_B_ISO(alpha_wt_2*180/np.pi, Possible_solutions_wolfrom[i][10], Possible_solutions_wolfrom[i][3], -Possible_solutions_wolfrom[i][0],m_t_2*10**(-3), x_P2, x_R2, Eps_alpha_2, 'int')
                        sigma_HP_R2 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP)/SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_wolfrom[i][3],-Possible_solutions_wolfrom[i][0], material, m_t_2*10**(-3), m_t_2*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)), alpha_wt_2)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[i][5])
                        sigma_HP_P2 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_wolfrom[i][3],-Possible_solutions_wolfrom[i][0], material, m_t_2*10**(-3), m_t_2*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)), alpha_wt_2)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[i][5])
                        sigma_HG_R2 = sigma_HP_R2*SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        sigma_HG_p2 = sigma_HP_P2*SH_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180))
                        if sigma_FG_R2/(sigma_R2*10**(-6)) > SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)) and  sigma_FG_p2/(sigma_p2*10**(-6)) > SF_min(m_t_2*np.cos(Possible_solutions_wolfrom[i][10]*np.pi/180)):
                            possibilities_2.append([(sigma_FG_R2/(sigma_R2*10**(-6))), sigma_FG_p2/(sigma_p2*10**(-6)), sigma_p2*10**(-6), sigma_R2*10**(-6), x_P2, x_R2, a_2,  sigma_H_R2*10**(-6), sigma_H_P2*10**(-6), sigma_HP_R2, sigma_HP_P2, sigma_HG_R2/(sigma_H_R2*10**(-6)), sigma_HG_p2/(sigma_H_P2*10**(-6))])
                            x_P2 += step
                            if x_P2 < underlimit:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
                        else:
                            x_P2 += step
                            if x_P2 < underlimit:
                                if len(possibilities_2) == 0:
                                    bending_stresses_after_PS[i].extend([0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
                                    exit_inner_loops_second_time = True
                                else:
                                    Best_config_2 = max(possibilities_2, key=lambda x:(abs(x[0]-x[1])))
                                    bending_stresses_after_PS[i].extend(Best_config_2)
                                    possibilities_2.clear()
                                    exit_inner_loops_second_time = True
                                break
            if exit_inner_loops_second_time:
                break
    return bending_stresses_after_PS
            
     
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################



#basic wolfrom calculations
def calculate_tangential_force_basic_wolfrom(Possible_solutions, values,mt1,mt2,NOP): # tangential forces are calculated by dividing the torque by the radius of the gear
    global Forces
    Forces = []
    Tr2 = float(values[1])
    for i in Possible_solutions:
        GR = i[4]
        eff = i[5]
        Fr2 = Tr2/((i[0]*mt2*10**-3)/2) #torque/radius
        Fp2 = Fr2/NOP
        Fr1 = (-(1-1/GR*1/eff)*Tr2)/((i[1]*mt1*10**-3)/2)
        Fp1 = Fr1/NOP
        Forces.append([round(Fr2,2), round(abs(Fr1),2), round(abs(Fp1),2), round(Fp2,2)])
    return Forces #in N

def calculate_bending_stress_basic_wolfrom(Forces, Possible_solutions_wolfrom, m_t_1, m_t_2,axial_length, NOP, HA_1, HA_2, Epsilon_alpha_2, Epsilon_alpha_1, type_of_gear, material, allowances):
    global bending_stress, permissible_bending_stress, S
    bending_stress = [] #for each gear the bending stress in this order: R2, R1,P1, P2
    permissible_bending_stress = []
    S = []
    for i in Possible_solutions_wolfrom:
        XE_1_ext = Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)) #generating profile shift coefficient for internal and external gears in each stage
        XE_2_ext = Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))    
        XE_1_int = Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        XE_2_int = Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        b1,b2 = calculate_width(Forces[Possible_solutions_wolfrom.index(i)][2], Forces[Possible_solutions_wolfrom.index(i)][3], axial_length, type_of_gear, m_t_1, m_t_2)
        Y_F_R2,rho_F_R2, sFn_R2, hFe_R2  = calculate_Y_F_ISO('int',i[0],HA_2[Possible_solutions_wolfrom.index(i)], XE_2_int, m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_2*10**-3*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180), Epsilon_alpha_2[Possible_solutions_wolfrom.index(i)], 0)
        Y_F_R1,rho_F_R1, sFn_R1, hFe_R1  = calculate_Y_F_ISO('int',i[1],HA_1[Possible_solutions_wolfrom.index(i)], XE_1_int, m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_1*10**-3*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180), Epsilon_alpha_1[Possible_solutions_wolfrom.index(i)], 0)
        Y_F_p1,rho_F_p1, sFn_p1, hFe_p1  = calculate_Y_F_ISO('ext',i[2],HA_1[Possible_solutions_wolfrom.index(i)], XE_1_ext, m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_1*10**-3*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180), Epsilon_alpha_1[Possible_solutions_wolfrom.index(i)], 0)
        Y_F_p2,rho_F_p2, sFn_p2, hFe_p2  = calculate_Y_F_ISO('ext',i[3],HA_2[Possible_solutions_wolfrom.index(i)], XE_2_ext, m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_2*10**-3*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180), Epsilon_alpha_2[Possible_solutions_wolfrom.index(i)], 0)
        sigma_R2 = Forces[Possible_solutions_wolfrom.index(i)][0]/NOP/(b2*10**-3*m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3)*Y_F_R2*calculate_Y_beta_ISO(overlap_ratio, HA_2[Possible_solutions_wolfrom.index(i)])*calculate_Y_B_ISO('int', m_t_2*10**-3*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180), 3.5*m_t_2*10**-3*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R2, rho_F_R2, hFe_R2)*K
        sigma_R1 = Forces[Possible_solutions_wolfrom.index(i)][1]/NOP/(b1*10**-3*m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3)*Y_F_R1*calculate_Y_beta_ISO(overlap_ratio, HA_1[Possible_solutions_wolfrom.index(i)])*calculate_Y_B_ISO('int', m_t_1*10**-3*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180), 3.5*m_t_1*10**-3*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R1, rho_F_R1, hFe_R1)*K
        sigma_p1 = Forces[Possible_solutions_wolfrom.index(i)][2]/(b1*10**-3*m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3)*Y_F_p1*calculate_Y_S_ISO(sFn_p1, rho_F_p1, hFe_p1)*calculate_Y_beta_ISO(overlap_ratio, HA_1[Possible_solutions_wolfrom.index(i)])*calculate_Y_B_ISO('ext', m_t_1*10**-3*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180), 3.5*m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3)*Y_DT*K
        sigma_p2 = Forces[Possible_solutions_wolfrom.index(i)][3]/(b2*10**-3*m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3)*Y_F_p2*calculate_Y_S_ISO(sFn_p2, rho_F_p2, hFe_p2)*calculate_Y_beta_ISO(overlap_ratio, HA_2[Possible_solutions_wolfrom.index(i)])*calculate_Y_B_ISO('ext', m_t_2*10**-3*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180), 3.5*m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)*10**-3)*Y_DT*K
        sigma_FP_R2 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP)*Y_delta_rel_T_ISO(sFn_R2/(2*rho_F_R2), Material_properties[material][5],calculate_Y_S_ISO(sFn_R2, rho_F_R2, hFe_R2),1)*Y_R_rel_T_iso(3e6*NOP, material)*Y_X/SF_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)) 
        sigma_FP_R1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP)*Y_delta_rel_T_ISO(sFn_R1/(2*rho_F_R1), Material_properties[material][5], calculate_Y_S_ISO(sFn_R1, rho_F_R1, hFe_R1),1)*Y_R_rel_T_iso(3e6*NOP, material)*Y_X/SF_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)) #Y_NT = 2.5 because of the static behavior of the ring gear 
        sigma_FP_P1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p1/(2*rho_F_p1), Material_properties[material][5], calculate_Y_S_ISO(sFn_p1, rho_F_p1, hFe_p1),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_FP_P2 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p2/(2*rho_F_p2), Material_properties[material][5], calculate_Y_S_ISO(sFn_p2, rho_F_p2, hFe_p2),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_FG_R2 = sigma_FP_R2*SF_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_FG_R1 = sigma_FP_R1*SF_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_FG_p1 = sigma_FP_P1*SF_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_FG_p2 = sigma_FP_P2*SF_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        bending_stress.append([sigma_R2*10**(-6), sigma_R1*10**(-6), sigma_p1*10**(-6), sigma_p2*10**(-6)]) #sigma_R2*10**(-6), sigma_R1*10**(-6), sigma_p1*10**(-6), sigma_p2*10**(-6)
        permissible_bending_stress.append([sigma_FP_R2, sigma_FP_R1, sigma_FP_P1, sigma_FP_P2]) # the permissible stresses for each gear
        S.append([sigma_FG_R2/(sigma_R2*10**(-6)), sigma_FG_R1/(sigma_R1*10**(-6)), sigma_FG_p1/(sigma_p1*10**(-6)), sigma_FG_p2/(sigma_p2*10**(-6))]) #safety factors for each gear
        B2.append(b2)
        B1.append(b1)
    return bending_stress, permissible_bending_stress,S #in MPa

def calculate_Hertzian_Stress_basic_wolfrom(Forces, Possible_solutions_wolfrom,m_t_1,m_t_2, NOP, HA_1, HA_2, alpha_n, material, x, B2, B1): #Hertzian stresses are calculated here 
    global Hertzian_stresses
    Hertzian_stresses = [] #for each gear the Hertzian stress in this order: first stage, second stage, third stage
    Permissible_Hertzian_stresses=[]
    safety_hertzian = []
    Z_E = calculate_ZE_ISO(Material_properties, material)
    for i in Possible_solutions_wolfrom:
        sigma_R2 = calculate_ZH_ISO(HA_2[Possible_solutions_wolfrom.index(i)], np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)), alpha_n, 'int')*Z_E*calculate_Z_eps_ISO(HA_2[Possible_solutions_wolfrom.index(i)], Epsilon_alpha_2[Possible_solutions_wolfrom.index(i)])*calculate_Z_beta_ISO(HA_2[Possible_solutions_wolfrom.index(i)])*np.sqrt((Forces[Possible_solutions_wolfrom.index(i)][0]/NOP)/(m_t_2*10**(-3)*i[3]*B2[Possible_solutions_wolfrom.index(i)]*10**(-3))*(-i[0]/i[3]+1)/(-i[0]/i[3]))*np.sqrt(K) #Z_D is = 1 for internal gears
        sigma_R1 = calculate_ZH_ISO(HA_1[Possible_solutions_wolfrom.index(i)], np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)), alpha_n, 'int')*Z_E*calculate_Z_eps_ISO(HA_1[Possible_solutions_wolfrom.index(i)], Epsilon_alpha_1[Possible_solutions_wolfrom.index(i)])*calculate_Z_beta_ISO(HA_1[Possible_solutions_wolfrom.index(i)])*np.sqrt((Forces[Possible_solutions_wolfrom.index(i)][1]/NOP)/(m_t_1*10**(-3)*i[2]*B1[Possible_solutions_wolfrom.index(i)]*10**(-3))*(-i[1]/i[2]+1)/(-i[1]/i[2]))*np.sqrt(K) #Z_D is = 1 for internal gears
        sigma_P1 = sigma_R1*calculate_Z_B_ISO(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))*180/np.pi, HA_1[Possible_solutions_wolfrom.index(i)], i[2], -i[1],m_t_1*10**(-3), x,0, Epsilon_alpha_1[Possible_solutions_wolfrom.index(i)], 'int')
        sigma_P2 = sigma_R2*calculate_Z_B_ISO(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))*180/np.pi, HA_2[Possible_solutions_wolfrom.index(i)], i[3], -i[0],m_t_2*10**(-3), x,0, Epsilon_alpha_2[Possible_solutions_wolfrom.index(i)], 'int')
        sigma_HP_R2 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP)/SH_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[3],-i[0], material, m_t_2*10**(-3), m_t_2*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[Possible_solutions_wolfrom.index(i)][5])
        sigma_HP_P2 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[3],-i[0], material, m_t_2*10**(-3), m_t_2*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[Possible_solutions_wolfrom.index(i)][5])
        sigma_HP_R1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP)/SH_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[2],-i[1], material, m_t_1*10**(-3), m_t_1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[Possible_solutions_wolfrom.index(i)][6])
        sigma_HP_P1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[2],-i[1], material, m_t_1*10**(-3), m_t_1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds[Possible_solutions_wolfrom.index(i)][6])
        sigma_HG_R2 = sigma_HP_R2*SH_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_HG_R1 = sigma_HP_R1*SH_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_HG_p1 = sigma_HP_P1*SH_min(m_t_1*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        sigma_HG_p2 = sigma_HP_P2*SH_min(m_t_2*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))
        Hertzian_stresses.append([sigma_R2*10**(-6), sigma_R1*10**(-6), sigma_P1*10**(-6), sigma_P2*10**(-6)])
        Permissible_Hertzian_stresses.append([sigma_HP_R2, sigma_HP_R1, sigma_HP_P1, sigma_HP_P2])
        safety_hertzian.append([sigma_HG_R2/(sigma_R2*10**(-6)), sigma_HG_R1/(sigma_R1*10**(-6)), sigma_HG_p1/(sigma_P1*10**(-6)), sigma_HG_p2/(sigma_P2*10**(-6))])
    return Hertzian_stresses, Permissible_Hertzian_stresses, safety_hertzian #in MPa

def calculate_speeds_basic_wolfrom(Possible_solutions, wR2, m1,m2):
    global speeds
    speeds = []
    for i in Possible_solutions:
        db_pl1 = m1*10**(-3)*i[2]*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)))
        db_pl2 = m2*10**(-3)*i[3]*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)))
        dbR1 = m1*10**(-3)*i[1]*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)))
        dbR2 = m2*10**(-3)*i[0]*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)))
        dapl_2 = m2*10**(-3)*i[3] + 2*m2*10**(-3)*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)
        dapl_1 = m1*10**(-3)*i[2]  + 2*m1*10**(-3)*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180) 
        daR_1 = m1*10**(-3)*i[1] - 2*m1*10**(-3)*np.cos(HA_1[Possible_solutions_wolfrom.index(i)]*np.pi/180)
        daR_2 = m2*10**(-3)*i[0] - 2*m2*10**(-3)*np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)

        dw_p2 = db_pl2/np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)))
        dw_ring2 = dbR2/np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180)))
        aw = 1/2*(dw_ring2 - dw_p2)
        dNf1_2 = np.sqrt((2*aw*np.sin(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))) - np.sqrt(daR_2**2 - dbR2**2))**2 + db_pl2**2) #planeet
        dNf2_2 = np.sqrt((2*aw*np.sin(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[Possible_solutions_wolfrom.index(i)]*np.pi/180))) + np.sqrt(dapl_2**2 - db_pl2**2))**2 + dbR2**2) #ring
        
        alpha_Ypl2 = np.arccos(db_pl2/dapl_2)
        alpha_Ypl1 = np.arccos(db_pl1/dapl_1)
        alpha_Y_R1 = np.arccos(dbR1/daR_1)
        alpha_Y_R2 = np.arccos(dbR2/daR_2)

        #rotational speeds in rpm
        wC =  wR2*i[4] 
        wR1 = -wC
        wP1 = (wR2-wR2*i[4])*i[0]/i[3]+wR2*i[4] - wC
        wP2 = (wR2-wR2*i[4])*i[0]/i[3]+wR2*i[4] - wC
        #tang velocities in m/s
        vt_R2 = wR2*np.pi/30*(m2*10**(-3)*i[0]/2)
        vt_R1 = wC*np.pi/30*(m1*10**(-3)*i[1]/2)
        vt_P1 = ((wR2*np.pi/30-wR2*np.pi/30*i[4])*i[0]/i[3] + wR2*np.pi/30*i[4])*(m1*10**(-3)*i[2]/2)
        vt_P2 = ((wR2*np.pi/30-wR2*np.pi/30*i[4])*i[0]/i[3]+ wR2*np.pi/30*i[4])*(m2*10**(-3)*i[3]/2)

        #tangential velocities in begin or end part of line of contact:
        VYt_pl1 = db_pl1/2 * wP1*np.pi/30 * np.tan(alpha_Ypl1)
        VYt_pl2 = db_pl2/2 * wP2*np.pi/30 * np.tan(alpha_Ypl2)
        VYt_R1 = dbR1/2 * wR1*np.pi/30 * np.tan(alpha_Y_R1)
        VYt_R2 = dbR2/2 * wR2*np.pi/30 * np.tan(alpha_Y_R2)
        sliding_1 = VYt_R1 - VYt_pl1
        sliding_2 = VYt_R2 - VYt_pl2
        speeds.append([wR2, wR1 , wC, wP1, wP2, vt_R2, vt_R1, vt_P1, vt_P2, VYt_R2, VYt_R1, VYt_pl1, VYt_pl2, sliding_1, sliding_2]) # rotational speeds: R2, R1, C, P1, P2 , tang speed: R2, R1, P1, P2, roll speeds: R2, R1, P1, P2
    return speeds #in rpm

def calculate_basic_wolfrom(): #main function that calculates all the possible solutions
    global values, Possible_solutions_wolfrom, overlap_ratio, sorted_possible_solutions_wolfrom, m_t_2, m_t_1, hole_decision, hole_diameter, pitch_diam
    Possible_solutions_wolfrom = []
    values = [entry.get() for entry in entry_boxes] #get a list of all the input values
    if module1_ok:
        m1 = float(values[7]) #is m_t in the case of helical gears
    else:
        if float(values[7]) < m_min_bending_stage_1*10**3:
            m1 = m_min_bending_stage_1*10**3
        else:
            m1 = float(values[7]) #is m_t in teh case of helical gears

    if module2_ok:
        m2 = float(values[8])#is m_t in the case of helical gears
    else:
        if float(values[8]) < m_min_bending_stage_2*10**3:
            m2 = m_min_bending_stage_2*10**3
        else:
            m2 = float(values[8]) 
   
    pitch_diam = int(values[0])
    GR_wolfrom = float(values[4])
    percentage_wolfrom = int(values[5])
    axial_length = float(values[3])
    Tr2 = float(values[1])
    material = material_var.get()
    Type_of_gear = gear_var.get()
    wolfrom_eff_setpoint = float(values[9])/100
    hole_decision = hole_var.get()
    if hole_decision =='Yes':
        hole_diameter = float(values[19])
    else:
        hole_diameter = 0
    allowances = allowance_var.get()
    wR2 = float(values[2])
    m_t_1 = m1
    m_t_2 = m2

    if percentage_wolfrom == "": #such that there is always a range of -10% and +10% if the user wouldn't enter anything.
        percentage_wolfrom = 10

    Zr1_max = int(pitch_diam/m_t_1) #maximum possible value
    Zr2_max = int(pitch_diam/m_t_2) 
    Number_of_planets = int(values[6])
    sin = np.sin(np.pi/Number_of_planets)

    #calcualting the wolfrom stage
    for Zr2 in range(Zr2_max, 22,-1):
        for Zr1 in range(Zr1_max, 22,-1):
            for Zp2 in range(int(Zr2 / 2 - hole_diameter/(2*m_t_2)),11,-1):
                for Zp1 in range(int(Zr1/2 - hole_diameter /(2*m_t_1)),11,-1):
                    S2 = Zr2/Zp2
                    S1 = Zr1/Zp1
                    if S1 != S2:
                        if  (GR_wolfrom-GR_wolfrom*percentage_wolfrom/100) <= (S2/(S2-S1)) <= (GR_wolfrom + GR_wolfrom*percentage_wolfrom/100): #chekc if gear ratio lies in the range
                            if Type_of_gear == 'Spur':
                                Helix_angle_input = 0           
                                Helix_angle_output = 0
                                overlap_ratio = 0
                                beta_b_1 = np.arctan(np.tan(Helix_angle_input*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input*np.pi/180))))
                                beta_b_2 = np.arctan(np.tan(Helix_angle_output*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output*np.pi/180))))
                                p_b_2_t = np.pi*m_t_2*10**(-3)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output*np.pi/180)))
                                p_b_1_t = np.pi*m_t_1*10**(-3)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input*np.pi/180)))
                                da2_2 = m_t_2*10**(-3)*Zr2 - 2*m_t_2*10**(-3) #for the internal one
                                db2_2 = m_t_2*10**(-3)*Zr2*np.cos(pressure_angle*math.pi/180) #for the internal one
                                da1_2 = m_t_2*10**(-3)*Zp2 + 2*m_t_2*10**(-3) #for the external one
                                db1_2 = m_t_2*10**(-3)*Zp2*np.cos(pressure_angle*math.pi/180) #for the external one
                                a_2 = (m_t_2*10**(-3)*Zr2 - m_t_2*10**(-3)*Zp2)/2 
                                #transverse contact ratio input mesh
                                da2_1 = m_t_1*10**(-3)*Zr1 - 2*m_t_1*10**(-3) 
                                db2_1= m_t_1*10**(-3)*Zr1*np.cos(pressure_angle*math.pi/180)
                                da1_1 = m_t_1*10**(-3)*Zp1  + 2*m_t_1*10**(-3) 
                                db1_1 = m_t_1*10**(-3)*Zp1*np.cos(pressure_angle*math.pi/180)
                                a_1 = (m_t_1*10**(-3)*Zr1 - m_t_1*10**(-3)*Zp1)/2  
                                Eps_alpha_2 = 1/p_b_2_t*(np.sqrt(da1_2**2 - db1_2**2)/2 - np.sqrt(da2_2**2 - db2_2**2)/2 + a_2*np.sin(pressure_angle*math.pi/180))
                                Eps_alpha_1 = 1/p_b_1_t*(np.sqrt(da1_1**2-db1_1**2)/2 - np.sqrt(da2_1**2-db2_1**2)/2 + a_1*np.sin(pressure_angle*math.pi/180))
                                Epsilon_approach_1 = -Zr1/(2*np.pi) * (np.tan(np.arccos(db2_1/da2_1)) - np.tan(pressure_angle*np.pi/180))
                                Epsilon_recess_1 = Zp1/(2*np.pi) * (np.tan(np.arccos(db1_1/da1_1)) - np.tan(pressure_angle*np.pi/180))
                                Epsilon_approach_2 = -Zr2/(2*np.pi) * (np.tan(np.arccos(db2_2/da2_2)) - np.tan(pressure_angle*np.pi/180))
                                Epsilon_recess_2 = Zp2/(2*np.pi) * (np.tan(np.arccos(db1_2/da1_2)) - np.tan(pressure_angle*np.pi/180))
                                f_R1_P1 = (fric_coeff*np.pi*(1/Zp1 - 1/Zr1)*(1-Eps_alpha_1 + Epsilon_approach_1**2 + Epsilon_recess_1**2))/np.cos(beta_b_1)
                                f_R2_P2 = (fric_coeff*np.pi*(1/Zp2 - 1/Zr2)*(1-Eps_alpha_2 + Epsilon_approach_2**2 + Epsilon_recess_2**2))/np.cos(beta_b_2) 
                                real_eff = (f_R1_P1+1)/(f_R1_P1*(S2/(S2-S1))+f_R2_P2*((S2/(S2-S1))-1)+1)
                            else:
                                if Type_of_gear == "Helical, eps_beta = 1":
                                    overlap_ratio = 1
                                if Type_of_gear == "Helical, eps_beta = 2":
                                    overlap_ratio = 2
                                error_1 = 5
                                error_2 = 5
                                Helix_angle_input = 0           
                                Helix_angle_output = 0
                                a_2 = (m_t_2*10**(-3)*Zr2 - m_t_2*10**(-3)*Zp2)/2
                                a_1 = (m_t_1*10**(-3)*Zr1 - m_t_1*10**(-3)*Zp1)/2  
                                beta_b_1_start = np.arctan(np.tan(Helix_angle_input*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input*np.pi/180))))
                                beta_b_2_start = np.arctan(np.tan(Helix_angle_output*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output*np.pi/180))))
                                while error_1 and error_2 > 10**(-8): #whiel loop that converges to a certain helix angle in order to have the correct overlap ratio
                                    da2_2 = m_t_2*10**(-3)*Zr2 - 2*m_t_2*10**(-3)*np.cos(Helix_angle_output)
                                    da1_2 = m_t_2*10**(-3)*Zp2 + 2*m_t_2*10**(-3)*np.cos(Helix_angle_output)
                                    da2_1 = m_t_1*10**(-3)*Zr1 - 2*m_t_1*10**(-3)*np.cos(Helix_angle_input) 
                                    da1_1 = m_t_1*10**(-3)*Zp1  + 2*m_t_1*10**(-3)*np.cos(Helix_angle_input) 
                                    db2_2 = m_t_2*10**(-3)*Zr2*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))) #for the internal one => must be alpha_t!
                                    db1_2 = m_t_2*10**(-3)*Zp2*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))) #for the external one => must be alpha_t!
                                    db2_1= m_t_1*10**(-3)*Zr1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))) #=> must be alpha_t!
                                    db1_1 = m_t_1*10**(-3)*Zp1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))) #=> must be alpha_t!
                                    Epsilon_approach_1 = -Zr1/(2*np.pi) * (np.tan(np.arccos(db2_1/da2_1)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))))
                                    Epsilon_recess_1 = Zp1/(2*np.pi) * (np.tan(np.arccos(db1_1/da1_1)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))))
                                    Epsilon_approach_2 = -Zr2/(2*np.pi) * (np.tan(np.arccos(db2_2/da2_2)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))))
                                    Epsilon_recess_2 = Zp2/(2*np.pi) * (np.tan(np.arccos(db1_2/da1_2)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))))
                                    p_b_2_t = np.pi*m_t_2*10**(-3)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output)))
                                    p_b_1_t = np.pi*m_t_1*10**(-3)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input)))
                                    Eps_alpha_2 = 1/p_b_2_t*(np.sqrt(da1_2**2 - db1_2**2)/2 - np.sqrt(da2_2**2 - db2_2**2)/2 + a_2*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))))
                                    Eps_alpha_1 = 1/p_b_1_t*(np.sqrt(da1_1**2-db1_1**2)/2 - np.sqrt(da2_1**2-db2_1**2)/2 + a_1*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))))
                                    f_R1_P1 = (fric_coeff*np.pi*(1/Zp1 - 1/Zr1)*(1-Eps_alpha_1 + Epsilon_approach_1**2 + Epsilon_recess_1**2))/np.cos(beta_b_1_start)
                                    f_R2_P2 = (fric_coeff*np.pi*(1/Zp2 - 1/Zr2)*(1-Eps_alpha_2 + Epsilon_approach_2**2 + Epsilon_recess_2**2))/np.cos(beta_b_2_start) 
                                    real_eff = (f_R1_P1+1)/(f_R1_P1*(S2/(S2-S1))+f_R2_P2*((S2/(S2-S1))-1)+1)
                                    Fr2 = Tr2/((Zr2*m_t_2*10**-3)/2) #torque/radius
                                    Ftp2 = Fr2/Number_of_planets  
                                    Fr1 = (-(1-1/(S2/(S2-S1))*1/real_eff)*Tr2)/((Zr1*m_t_1*10**-3)/2)
                                    Ftp1 = Fr1/Number_of_planets
                                    b1, b2 = calculate_width(Ftp1,Ftp2, axial_length, Type_of_gear, m_t_1, m_t_2)
                                    Helix_angle_input = np.arctan(overlap_ratio*np.pi*m_t_1/b1)          
                                    Helix_angle_output = np.arctan(overlap_ratio*np.pi*m_t_2/b2)
                                    beta_b_1 = np.arctan(np.tan(Helix_angle_input)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input))))
                                    beta_b_2 = np.arctan(np.tan(Helix_angle_output)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output))))
                                    error_1 = abs(beta_b_1 - beta_b_1_start)
                                    error_2 = abs(beta_b_2 - beta_b_2_start)
                                    beta_b_1_start = beta_b_1
                                    beta_b_2_start = beta_b_2
                            fR1_P1_ = fric_coeff*np.pi*(1/Zp1 - 1/Zr1)*0.5*np.cos(beta_b_1)
                            fR2_P2_ = fric_coeff*np.pi*(1/Zp2 - 1/Zr2)*0.5/np.cos(beta_b_2)
                            eff = (fR1_P1_+1)/(fR1_P1_*(S2/(S2-S1))+fR2_P2_*((S2/(S2-S1))-1)+1)
                            if  eff>= wolfrom_eff_setpoint: #alle bij houden die goed genoeg zijn met eps_alplha = 1, en de rest 0.5
                                if abs((Zr1*m_t_1 - Zp1*m_t_1) - (Zr2*m_t_2 -Zp2*m_t_2)) > min(m_t_1, m_t_2): #coaxiality 
                                    continue
                                if hole_decision == 'Yes':
                                    if (Zr2*m_t_2/2-2*(Zp2*m_t_2/2) - m_t_2*np.cos(Helix_angle_output)) < hole_diameter/2 + margin: #clearance/2 want je moet maar langs 1 kant clearance hebben
                                        continue
                                    if (Zr1*m_t_1/2 -2*(Zp1*m_t_1/2) - m_t_1*np.cos(Helix_angle_input)) < hole_diameter/2 + margin: #neighbouring voor cable hole
                                        continue    
                                if (Zr1 - Zp1)*sin < Zp1 + 2*np.cos(Helix_angle_input) + clearance/(m_t_1): #neighbouring
                                    continue
                                if (Zr2 - Zp2)*sin < Zp2 + 2*np.cos(Helix_angle_output) + clearance/(m_t_2): #neighbouring
                                    continue
                                if (Zp1*Zr2-Zp2*Zr1)/Number_of_planets != int((Zp1*Zr2-Zp2*Zr1)/Number_of_planets): #assembly condition, when identical planets are assumed
                                    continue
                                else:
                                    Possible_solutions_wolfrom.append([Zr2, Zr1, Zp1, Zp2, (S2/(S2-S1)), real_eff, eff]) #store every value if it was a good solution
                                    Epsilon_alpha_2.append(Eps_alpha_2)
                                    Epsilon_alpha_1.append(Eps_alpha_1)
                                    HA_1.append(Helix_angle_input*180/np.pi)
                                    HA_2.append(Helix_angle_output*180/np.pi)
                                    A.append([a_1,a_2])
                                    f_R2.append(f_R2_P2)
                                    f_R1.append(f_R1_P1)

    Forces = calculate_tangential_force_basic_wolfrom(Possible_solutions_wolfrom, values,m_t_1, m_t_2, Number_of_planets) #for all the solutions, calculate tangential forces, bending stresses, safety factors,...
    bending_stresses, permissible_bending_stresses, SF = calculate_bending_stress_basic_wolfrom(Forces, Possible_solutions_wolfrom, m_t_1, m_t_2,axial_length, Number_of_planets, HA_1, HA_2, Epsilon_alpha_2, Epsilon_alpha_1, Type_of_gear, material, allowances)
    rotational_velocity = calculate_speeds_basic_wolfrom(Possible_solutions_wolfrom, wR2,m_t_1,m_t_2)
    Hertzian_stresses, permissible_H_stresses,safety_H = calculate_Hertzian_Stress_basic_wolfrom(Forces, Possible_solutions_wolfrom,m_t_1,m_t_2, Number_of_planets, HA_1, HA_2,pressure_angle, material,0, B2, B1)

    for i in range(len(Forces)): #append all the resutlts to each solution
        Possible_solutions_wolfrom[i].append(B1[i])
        Possible_solutions_wolfrom[i].append(B2[i])
        Possible_solutions_wolfrom[i].append(HA_1[i])
        Possible_solutions_wolfrom[i].append(HA_2[i])
        Possible_solutions_wolfrom[i].append(Epsilon_alpha_1[i])
        Possible_solutions_wolfrom[i].append(Epsilon_alpha_2[i])
        Possible_solutions_wolfrom[i].append(A[i][0])
        Possible_solutions_wolfrom[i].append(A[i][1])
        Possible_solutions_wolfrom[i].append(1-f_R1[i])
        Possible_solutions_wolfrom[i].append(1-f_R2[i])
        Possible_solutions_wolfrom[i].append(20)
        Possible_solutions_wolfrom[i].append(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[i]*np.pi/180))*180/np.pi)
        Possible_solutions_wolfrom[i].append(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[i]*np.pi/180))*180/np.pi)
        Possible_solutions_wolfrom[i].append(Forces[i][0])
        Possible_solutions_wolfrom[i].append(Forces[i][1])
        Possible_solutions_wolfrom[i].append(Forces[i][2])
        Possible_solutions_wolfrom[i].append(Forces[i][3])
        Possible_solutions_wolfrom[i].append(bending_stresses[i][0])
        Possible_solutions_wolfrom[i].append(bending_stresses[i][1])
        Possible_solutions_wolfrom[i].append(bending_stresses[i][2])
        Possible_solutions_wolfrom[i].append(bending_stresses[i][3])
        Possible_solutions_wolfrom[i].append(Hertzian_stresses[i][0])
        Possible_solutions_wolfrom[i].append(Hertzian_stresses[i][1])
        Possible_solutions_wolfrom[i].append(Hertzian_stresses[i][2])
        Possible_solutions_wolfrom[i].append(Hertzian_stresses[i][3])
        Possible_solutions_wolfrom[i].append(permissible_bending_stresses[i][0])
        Possible_solutions_wolfrom[i].append(permissible_bending_stresses[i][1])
        Possible_solutions_wolfrom[i].append(permissible_bending_stresses[i][2])
        Possible_solutions_wolfrom[i].append(permissible_bending_stresses[i][3])
        Possible_solutions_wolfrom[i].append(SF[i][0])
        Possible_solutions_wolfrom[i].append(SF[i][1])
        Possible_solutions_wolfrom[i].append(SF[i][2])
        Possible_solutions_wolfrom[i].append(SF[i][3])
        Possible_solutions_wolfrom[i].append(permissible_H_stresses[i][0])
        Possible_solutions_wolfrom[i].append(permissible_H_stresses[i][1])
        Possible_solutions_wolfrom[i].append(permissible_H_stresses[i][2])
        Possible_solutions_wolfrom[i].append(permissible_H_stresses[i][3])
        Possible_solutions_wolfrom[i].append(safety_H[i][0])
        Possible_solutions_wolfrom[i].append(safety_H[i][1])
        Possible_solutions_wolfrom[i].append(safety_H[i][2])
        Possible_solutions_wolfrom[i].append(safety_H[i][3]) #47
  
    sorted_possible_solutions_wolfrom = sorted(Possible_solutions_wolfrom, key=lambda x: x[6], reverse = True) #sort the solution in terms of efficiency with transverse contact ratio =1
    sorted_possible_solutions_wolfrom = sorted_possible_solutions_wolfrom[:20] #keep only the 20 best solutions to perform a profile shift
    if len(Possible_solutions_wolfrom) < 20: #if len(list) smaller than 20, use all the solutions
        sorted_possible_solutions_wolfrom = sorted_possible_solutions_wolfrom[:len(Possible_solutions_wolfrom)]
    Na_PS = profile_shift_basic_wolfrom(material, sorted_possible_solutions_wolfrom, m_t_1, m_t_2, Number_of_planets, allowances) #perform profile shift on all these 20 or less solutions
 
    for i in range(len(sorted_possible_solutions_wolfrom)): #store all the results of the profile shift operation for each solution
        for j in range(24):
            sorted_possible_solutions_wolfrom[i].append(Na_PS[i][j])

    for i in range(len(sorted_possible_solutions_wolfrom)):  
        if sorted_possible_solutions_wolfrom[i][65] > pitch_diam*10**(-3)/2 - sorted_possible_solutions_wolfrom[i][3]*m_t_2*10**(-3)/2:
            sorted_possible_solutions_wolfrom[i].append('After profile shift, the solution does not fit anymore, this is not a good solution!')
        else:
            sorted_possible_solutions_wolfrom[i].append('After profile shift, the solution still fits!')
        
        if sorted_possible_solutions_wolfrom[i][65] == 0 or sorted_possible_solutions_wolfrom[i][54] == 0:
            sorted_possible_solutions_wolfrom[i].append('This is not a good solution, the profile shift could not be used to comply with the coaxiality condition, to prevent undercutting or to make the gears strong enough.')


    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, f"There are {len(Possible_solutions_wolfrom)} solutions", "bold") #give a message of how many solutions there are
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "Possible configurations for basic Wolfrom part: ", "bold")
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "{:<7}{:<5}{:<5}{:<5}{:<5}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:15}\n".format("#", "Zr2", "Zr1", "Zp1", "Zp2", "Wolfrom_GR", "Wolfrom_eff","eff_ea=1", "b1","b2","HA_1","HA_2","Eps_a_1", "Eps_a_2", 'a_1', 'a_2',
                                                                                                                                                                                                                                                               'mesh_eff_R1', 'mesh_eff_R2', 'alpha_n', 'alpha_t_1', 'alpha_t_2', 'Ft_R2', 'Ft_R1', 'Ft_p1', 'Ft_p2', 'BS_R2', 'BS_R1', 'BS_p1', 'BS_p2', 'HS_R2', 'HS_R1', 
                                                                                                                                                                                                                                                               'HS_p1', 'HS_p2', 'BS_perm_R2', 'BS_perm_R1', 'BS_perm_p1', 'BS_perm_p2', 'SF_R2', 'SF_R1', 'SF_p1', 'SF_p2', 'H_R2_perm', 'H_R1_perm', 'H_P1_perm', 'H_P2_perm', 'SH_R2', 'SH_R1', 'SH_P1', 'SH_P2','SF_x_R1', 'SF_x_P1', 'BS_P1', 'BS_R1', 'x_P1', 'x_R1', 'a_1_na_PS','sigma_H_R1',' sigma_H_P1', 'S_H_x_R1', 'S_H_x_P1', 'SF_x_R2', 'SF_x_P2', 'BS_P2', 'BS_R2', 'x_P2', 'x_R2', 'a_2_After_PS', 'sigma_H_R2', 'sigma_H_P2', 'sigma_perm_HR2', 'sigma_perm_HP2', 'S_H_x_R2', 'S_H_x_P2', 'Warning'))

    # Loop through Possible_solutions_wolfrom and display the values in the interface
    for solution in sorted_possible_solutions_wolfrom[:]: 
        output_text.insert(tk.END, "{:<7}{:<5}{:<5}{:<5}{:<5}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{}\n".format(
        sorted_possible_solutions_wolfrom.index(solution) + 1, solution[0], solution[1], solution[2], solution[3], round(solution[4], 3), 
        round(solution[5], 5), round(solution[6], 4), round(solution[7], 4), round(solution[8], 4), round(solution[9], 4), 
        round(solution[10], 4), round(solution[11], 4), round(solution[12], 3), round(solution[13], 3), round(solution[14], 3), 
        round(solution[15], 3), round(solution[16], 3), round(solution[17], 3), round(solution[18], 3), round(solution[19], 3), round(solution[20], 3), round(solution[21], 3), round(solution[22], 3), 
        round(solution[23], 3), round(solution[24], 3), round(solution[25], 3) ,round(solution[26], 3), round(solution[27], 3),round(solution[28], 3),round(solution[29], 3),round(solution[30], 3),
        round(solution[31], 3),round(solution[32], 3),round(solution[33], 3),round(solution[34], 3), round(solution[35], 3), round(solution[36], 3), round(solution[37], 3), round(solution[38], 3), round(solution[39], 3), 
        round(solution[40],3) , round(solution[41],3) , round(solution[42],3) , round(solution[43],3), round(solution[44], 5), round(solution[45], 5), round(solution[46], 5), round(solution[47], 5), round(solution[48], 5),
        round(solution[49], 5),round(solution[50], 5),round(solution[51], 5), round(solution[52],3) , round(solution[53],5) , round(solution[54],5) , round(solution[55],5) , round(solution[56], 5), round(solution[57], 5), 
        round(solution[58], 5), round(solution[59], 5), round(solution[60], 5), round(solution[61], 5), round(solution[62], 5), round(solution[63], 5), round(solution[64], 5), round(solution[65], 5), round(solution[66], 5), round(solution[67], 5), round(solution[68], 5), round(solution[69], 5), round(solution[70], 5), round(solution[71], 5), solution[72]))
    
    # Ensure the Text widget scrolls to the bottom to show new output
    output_text.yview(tk.END)
    output_text.insert(tk.END, "{:<7}{:<5}{:<5}{:<5}{:<5}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}{:<15}\n".format("#", "Zr2", "Zr1", "Zp1", "Zp2", "Wolfrom_GR", "Wolfrom_eff","eff_ea=1", "b1","b2","HA_1","HA_2","Eps_a_1", "Eps_a_2", 'a_1', 'a_2',
                                                                                                                                                                                                                                                               'mesh_eff_R1', 'mesh_eff_R2', 'alpha_n', 'alpha_t_1', 'alpha_t_2', 'Ft_R2', 'Ft_R1', 'Ft_p1', 'Ft_p2', 'BS_R2', 'BS_R1', 'BS_p1', 'BS_p2', 'HS_R2', 'HS_R1', 
                                                                                                                                                                                                                                                               'HS_p1', 'HS_p2', 'BS_perm_R2', 'BS_perm_R1', 'BS_perm_p1', 'BS_perm_p2', 'SF_R2', 'SF_R1', 'SF_p1', 'SF_p2', 'H_R2_perm', 'H_R1_perm', 'H_P1_perm', 'H_P2_perm', 'SH_R2', 'SH_R1', 'SH_P1', 'SH_P2','SF_x_R1', 'SF_x_P1', 'BS_P1', 'BS_R1', 'x_P1', 'x_R1', 'a_1_na_PS','sigma_H_R1',' sigma_H_P1', 'S_H_x_R1', 'S_H_x_P1', 'SF_x_R2', 'SF_x_P2', 'BS_P2', 'BS_R2', 'x_P2', 'x_R2', 'a_2_After_PS', 'sigma_H_R2', 'sigma_H_P2', 'sigma_perm_HR2', 'sigma_perm_HP2', 'S_H_x_R2', 'S_H_x_P2'))

    output_text.insert(tk.END, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", "bold")
    output_text.insert(tk.END, "Pick a solution and fill in the configuration number (the number in the first column) in the newly shown input field. Afterwards, click on the 'Calculate_pregearing' button to start the calculation for the pregearing.", "bold")
    show_configuration_number()

################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################
################################################################################################################################################################################################

#ALL THE SAME FUNCTIONS BUT HERE FOR THE PREGEARING



def calculate_tangential_force_pregearing(Possible_solutions, values, m_t_0, m_t__1, wolfrom_configuration_chosen, NOP_pre): # tangential forces are calculated by dividing the torque by the radius of the gear
    Forces = []
    GR_w = wolfrom_configuration_chosen[4]
    eff_w= wolfrom_configuration_chosen[5]
    for i in Possible_solutions:
        GR_0 = i[4]
        eff_tot = i[5]*eff_w
        Tr2 = float(values[1]) #given by the user
        Fs = (-1/eff_tot * 1/(GR_0*GR_w)*Tr2)/(m_t_0*10**-3*i[2]/2)
        Fp0 = Fs/NOP_pre
        Fr_1 = ((1/(GR_0*GR_w)*1/eff_tot + 1/GR_w*1/eff_w)*Tr2)/(m_t__1*10**-3*i[0]/2)
        Fp_1 = Fr_1/NOP_pre
        Forces.append([abs(Fr_1), abs(Fp_1), abs(Fs), abs(Fp0)])
    return Forces #in N

def calculate_speeds_pregearing(Possible_solutions, wolfrom_configuration_chosen, wR2, m_t_0, m_t__1):
    speeds= []
    wR2 = wR2*np.pi/30                                                                                                                                                                                                                                                                                                                                                                                                                                                                         
    for i in Possible_solutions:
        speeds.append([0,wR2*wolfrom_configuration_chosen[4],wolfrom_configuration_chosen[4]*i[4]*wR2,(wR2-wR2*wolfrom_configuration_chosen[4])*wolfrom_configuration_chosen[0]/wolfrom_configuration_chosen[3]+wR2*wolfrom_configuration_chosen[4], (wR2-wR2*wolfrom_configuration_chosen[4])*wolfrom_configuration_chosen[0]/wolfrom_configuration_chosen[3]+wR2*wolfrom_configuration_chosen[4], wolfrom_configuration_chosen[4]*i[4]*wR2*(i[2]*m_t_0*10**(-3)/2), wolfrom_configuration_chosen[4]*wR2*(i[0]*m_t__1*10**(-3)/2) ]) #wR_1
        # speeds.extend([wR2*wolfrom_configuration_chosen[4]]) #wc
        # speeds.extend([wolfrom_configuration_chosen[4]*i[4]*wR2]) #ws
        # speeds.extend([(wR2-wR2*wolfrom_configuration_chosen[4])*wolfrom_configuration_chosen[0]/wolfrom_configuration_chosen[3]+wR2*wolfrom_configuration_chosen[4]]) #wp_1
        # speeds.extend([(wR2-wR2*wolfrom_configuration_chosen[4])*wolfrom_configuration_chosen[0]/wolfrom_configuration_chosen[3]+wR2*wolfrom_configuration_chosen[4]]) #wp0

        # speeds.extend([wolfrom_configuration_chosen[4]*i[4]*wR2*(i[2]*m_t_0*10**(-3)/2)]) #v_t_s
    return speeds #in rad/s

def calculate_bending_stress_pregearing(Possible_solutions, m_t_0, m_t__1, Forces, NOP_pre, HA_0, HA__1, axial_length, Type_of_gear, material, allowances):
    bending_stress = [] 
    permissible_bending_stress = []
    SF = []
    for i in Possible_solutions:
        XE_0_ext = Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))
        XE__1_ext = Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        XE__1_int = Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))        
        b0, b_1 = calculate_width(Forces[Possible_solutions.index(i)][3], Forces[Possible_solutions.index(i)][1], axial_length, Type_of_gear, m_t_0, m_t__1)
        Y_F_R_1,rho_F_R_1, sFn_R_1, hFe_R_1  = calculate_Y_F_ISO('int',i[0], HA__1[Possible_solutions.index(i)], XE__1_int, m_t__1*10**-3*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180), pressure_angle,0.38*m_t__1*10**-3*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180), Epsilon_alpha__1[Possible_solutions.index(i)],0)
        Y_F_P_1,rho_F_P_1, sFn_P_1, hFe_P_1  = calculate_Y_F_ISO('ext',i[1], HA__1[Possible_solutions.index(i)], XE__1_ext, m_t__1*10**-3*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180), pressure_angle,0.38*m_t__1*10**-3*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180), Epsilon_alpha__1[Possible_solutions.index(i)],0)
        Y_F_S,rho_F_S, sFn_S, hFe_S  = calculate_Y_F_ISO('ext',i[2],HA_0[Possible_solutions.index(i)], XE_0_ext, m_t_0*10**-3*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180), pressure_angle,0.38*m_t_0*10**-3*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180), Epsilon_alpha_0[Possible_solutions.index(i)],0)
        Y_F_p0,rho_F_p0, sFn_p0, hFe_p0  = calculate_Y_F_ISO('ext',i[3],HA_0[Possible_solutions.index(i)], XE_0_ext, m_t_0*10**-3*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180), pressure_angle,0.38*m_t_0*10**-3*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180), Epsilon_alpha_0[Possible_solutions.index(i)],0)
        sigma_R_1 = Forces[Possible_solutions.index(i)][0]/NOP_pre/(b_1*10**(-3)*m_t__1*10**(-3)*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*Y_F_R_1*calculate_Y_S_ISO(sFn_R_1, rho_F_R_1, hFe_R_1)*calculate_Y_beta_ISO(overlap_ratio, HA__1[Possible_solutions.index(i)])*calculate_Y_B_ISO('int', m_t__1*10**-3*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180), 3.5*m_t__1*10**-3*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*Y_DT*K
        sigma_P_1 = Forces[Possible_solutions.index(i)][1]/(b_1*10**(-3)*m_t__1*10**(-3)*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*Y_F_P_1*calculate_Y_S_ISO(sFn_P_1, rho_F_P_1, hFe_P_1)*calculate_Y_beta_ISO(overlap_ratio,HA__1[Possible_solutions.index(i)])*calculate_Y_B_ISO('ext', m_t__1*10**(-3)*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180), 3.5*m_t__1*10**(-3)*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*Y_DT*K
        sigma_s = Forces[Possible_solutions.index(i)][2]/NOP_pre/(b0*10**(-3)*m_t_0*10**(-3)*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*Y_F_S*calculate_Y_S_ISO(sFn_S, rho_F_S, hFe_S)*calculate_Y_beta_ISO(overlap_ratio, HA_0[Possible_solutions.index(i)])*calculate_Y_B_ISO('ext', m_t_0*10**(-3)*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180), 3.5*m_t_0*10**(-3)*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*Y_DT*K
        sigma_P0 = Forces[Possible_solutions.index(i)][3]/(b0*10**(-3)*m_t_0*10**(-3)*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*Y_F_p0*calculate_Y_S_ISO(sFn_p0, rho_F_p0, hFe_p0)*calculate_Y_beta_ISO(overlap_ratio, HA_0[Possible_solutions.index(i)])*calculate_Y_B_ISO('ext', m_t_0*10**(-3)*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180), 3.5*m_t_0*10**(-3)*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*Y_DT*K
        sigma_FP_R_1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP_pre)*Y_delta_rel_T_ISO(sFn_R_1/(2*rho_F_R_1), Material_properties[material][5], calculate_Y_S_ISO(sFn_R_1, rho_F_R_1, hFe_R_1),1)*Y_R_rel_T_iso(3e6*NOP_pre, material)*Y_X/SF_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        sigma_FP_P_1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_P_1/(2*rho_F_P_1), Material_properties[material][5], calculate_Y_S_ISO(sFn_P_1, rho_F_P_1, hFe_P_1),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        sigma_FP_s = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP_pre)*Y_delta_rel_T_ISO(sFn_S/(2*rho_F_S), Material_properties[material][5], calculate_Y_S_ISO(sFn_S, rho_F_S, hFe_S),1)*Y_R_rel_T_iso(3e6*NOP_pre, material)*Y_X/SF_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))
        sigma_FP_P0 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p0/(2*rho_F_p0), Material_properties[material][5], calculate_Y_S_ISO(sFn_p0, rho_F_p0, hFe_p0),1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))
        sigma_FG_R_1 = sigma_FP_R_1*SF_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        sigma_FG_P_1 = sigma_FP_P_1*SF_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        sigma_FG_s = sigma_FP_s*SF_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))
        sigma_FG_P0 = sigma_FP_P0*SF_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))
        bending_stress.append([sigma_R_1*10**(-6), sigma_P_1*10**(-6), sigma_s*10**(-6), sigma_P0*10**(-6)])
        permissible_bending_stress.append([sigma_FP_R_1, sigma_FP_P_1, sigma_FP_s, sigma_FP_P0])
        SF.append([sigma_FG_R_1/(sigma_R_1*10**(-6)), sigma_FG_P_1/(sigma_P_1*10**(-6)), sigma_FG_s/(sigma_s*10**(-6)), sigma_FG_P0/(sigma_P0*10**(-6))])
        B0.append(b0)
        B_1.append(b_1)
    return bending_stress, permissible_bending_stress, SF #in MPa


def calculate_Hertzian_Stress_pregearing(Possible_solutions, m_t_0, m_t__1, Forces, NOP_pre, HA_0, HA__1, alpha_n, material, x, B0, B_1, speeds_pre): #Hertzian stresses are calculated here 
    Hertzian_stresses = [] 
    permissible_H_pregearing = []
    safety_H_pre = []
    Z_E = calculate_ZE_ISO(Material_properties, material)
    for i in Possible_solutions:
        sigma_R_1 = calculate_ZH_ISO(HA__1[Possible_solutions.index(i)], np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180)), alpha_n, 'int')*Z_E*calculate_Z_eps_ISO(HA__1[Possible_solutions.index(i)], Epsilon_alpha__1[Possible_solutions.index(i)])*calculate_Z_beta_ISO(HA__1[Possible_solutions.index(i)])*np.sqrt((Forces[Possible_solutions.index(i)][0]/NOP_pre)/(m_t__1*10**(-3)*i[1]*B_1[Possible_solutions.index(i)]*10**(-3))*(-i[0]/i[1]+1)/(-i[0]/i[1]))*np.sqrt(K) #Z_D is = 1 for internal gears
        sigma_s = calculate_ZH_ISO(HA_0[Possible_solutions.index(i)], np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)), alpha_n, 'ext')*Z_E*calculate_Z_eps_ISO(HA_0[Possible_solutions.index(i)], Epsilon_alpha_0[Possible_solutions.index(i)])*calculate_Z_beta_ISO(HA_0[Possible_solutions.index(i)])*np.sqrt((Forces[Possible_solutions.index(i)][2]/NOP_pre)/(m_t_0*10**(-3)*i[2]*B0[Possible_solutions.index(i)]*10**(-3))*(i[3]/i[2]+1)/(i[3]/i[2]))*np.sqrt(K)*calculate_Z_B_ISO(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*180/np.pi, HA_0[Possible_solutions.index(i)], i[2], i[3], m_t_0*10**(-3), x,0, Epsilon_alpha_0[Possible_solutions.index(i)], 'ext')
        sigma_P_1 = sigma_R_1*calculate_Z_B_ISO(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*180/np.pi, HA__1[Possible_solutions.index(i)], i[1], -i[0], m_t__1*10**(-3), 0,0, Epsilon_alpha__1[Possible_solutions.index(i)], 'int')
        sigma_P0 = calculate_ZH_ISO(HA_0[Possible_solutions.index(i)], np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)), alpha_n, 'ext')*Z_E*calculate_Z_eps_ISO(HA_0[Possible_solutions.index(i)], Epsilon_alpha_0[Possible_solutions.index(i)])*calculate_Z_beta_ISO(HA_0[Possible_solutions.index(i)])*np.sqrt((Forces[Possible_solutions.index(i)][3])/(m_t_0*10**(-3)*i[2]*B0[Possible_solutions.index(i)]*10**(-3))*(i[3]/i[2]+1)/(i[3]/i[2]))*np.sqrt(K)*calculate_Z_D_ISO(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)), HA_0[Possible_solutions.index(i)], i[2], i[3],m_t_0*10**(-3), x, Epsilon_alpha_0[Possible_solutions.index(i)])  
      
        sigma_HP_R_1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP_pre)/SH_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[1],-i[0], material, m_t__1*10**(-3), m_t__1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[Possible_solutions.index(i)][6])
        sigma_HP_P_1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[1],-i[0], material, m_t__1*10**(-3), m_t__1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[Possible_solutions.index(i)][6])
        sigma_HP_s = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP_pre)/SH_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[3],i[2], material, m_t_0*10**(-3), m_t_0*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[Possible_solutions.index(i)][5])
        sigma_HP_P0 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))*calculate_Z_R_ISO(i[3],i[2], material, m_t_0*10**(-3), m_t_0*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180)))*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[Possible_solutions.index(i)][5])
       
        sigma_HG_R_1 = sigma_HP_R_1*SH_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        sigma_HG_P_1 = sigma_HP_P_1*SH_min(m_t__1*np.cos(HA__1[Possible_solutions.index(i)]*np.pi/180))
        sigma_HG_s = sigma_HP_s*SH_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))
        sigma_HG_p0 = sigma_HP_P0*SH_min(m_t_0*np.cos(HA_0[Possible_solutions.index(i)]*np.pi/180))

        Hertzian_stresses.append([sigma_R_1*10**(-6), sigma_P_1*10**(-6), sigma_s*10**(-6), sigma_P0*10**(-6)])   
        permissible_H_pregearing.append([sigma_HP_R_1, sigma_HP_P_1, sigma_HP_s, sigma_HP_P0])
        safety_H_pre.append([sigma_HG_R_1/(sigma_R_1*10**(-6)), sigma_HG_P_1/(sigma_P_1*10**(-6)), sigma_HG_s/(sigma_s*10**(-6)), sigma_HG_p0/(sigma_P0*10**(-6))])     
    return Hertzian_stresses, permissible_H_pregearing, safety_H_pre


def profile_shift_pregearing(material, Possible_solutions_pregearing, m_t_0, m_t__1, NOP_pre, allowances, speeds_pre):
    possibilities__1 = []
    possibilities_0 = []
    bending_stresses_after_PS_0 = []
    Z_E = calculate_ZE_ISO(Material_properties, material)
    for i in range(len(Possible_solutions_pregearing)):
        exit_inner_loops = False
        exit_inner_loops_second_time = False
        db2__1 = m_t__1*10**(-3)*Possible_solutions_pregearing[i][0]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) # ring gear in output stage
        db1__1 = m_t__1*10**(-3)*Possible_solutions_pregearing[i][1]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) # planet in output stage
        db2_0 = m_t_0*10**(-3)*Possible_solutions_pregearing[i][2]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))) # sun in input stage
        db1_0 = m_t_0*10**(-3)*Possible_solutions_pregearing[i][3]*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))) # planet in input stage
        tau__1 = (2*np.pi)/(Possible_solutions_pregearing[i][1]) #for output stage
        tau_0 = (2*np.pi)/(Possible_solutions_pregearing[i][3]) #for input stage
        diameters = [wolfrom_configuration_chosen[0]*m_t_2, wolfrom_configuration_chosen[1]*m_t_1, Possible_solutions_pregearing[i][0]*m_t__1]
        max_diam = max(diameters)
        max_index = diameters.index(max_diam)
        if max_index == 0:
            m_max = m_t_2
            beta_max = wolfrom_configuration_chosen[10]
        elif max_index == 1:
            m_max = m_t_1
            beta_max = wolfrom_configuration_chosen[9]
        else:
            m_max = m_t__1
            beta_max = Possible_solutions_pregearing[i][10]
        y_max_0 = (max_diam/2*10**(-3))/(m_t_0*10**(-3) * np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) + (1.25*m_max*10**(-3)*np.cos(beta_max*np.pi/180))/(m_t_0*10**(-3)*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) - 1 + (3.5*m_max*10**(-3))/(m_t_0*10**(-3)*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) - (Possible_solutions_pregearing[i][2] + Possible_solutions_pregearing[i][3])/(2*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) - (Possible_solutions_pregearing[i][3])/(2*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
        cos_alphawt_0 = np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)))/((2*y_max_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)/(Possible_solutions_pregearing[i][2]+Possible_solutions_pregearing[i][3])) + 1)
        alpha_wt_0 = np.arccos(cos_alphawt_0)
        inv_alpha_wt_min_0 = np.tan(alpha_wt_0) - alpha_wt_0
        max_sum_of_prof_shift_0 = (inv_alpha_wt_min_0-np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))) + np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)))/(2*np.tan(pressure_angle*np.pi/180))*(Possible_solutions_pregearing[i][2] + Possible_solutions_pregearing[i][3])
        x_s = 1-(Possible_solutions_pregearing[i][2]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
        if x_s <0:
            x_s = 0
        while x_s <= max_sum_of_prof_shift_0:
            x_P0 = max_sum_of_prof_shift_0 - x_s
            if x_P0 < 1-(Possible_solutions_pregearing[i][3]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)):
                x_s +=0.005
                if x_s > max_sum_of_prof_shift_0:
                    if len(possibilities_0)==0:
                        bending_stresses_after_PS_0.append([0,0,0,0,0,0,0,0,0,0,0])
                        exit_inner_loops = True
                    else:
                        Best_config_0 = max(possibilities_0, key=lambda x:(x[0]))
                        bending_stresses_after_PS_0.append(Best_config_0)
                        possibilities_0.clear()
                        exit_inner_loops = True
                        exit_inner_loops_second_time = False
                    break
                continue
            XE_P0_ext = x_P0 + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            XE_s_ext = x_s + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            da2_0 = m_t_0*10**(-3)*Possible_solutions_pregearing[i][2] + 2*m_t_0*10**(-3)*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*(1+x_s)
            da1_0 = m_t_0*10**(-3)*Possible_solutions_pregearing[i][3] + 2*m_t_0*10**(-3)*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*(1+x_P0)
            dw_p0 = db1_0/np.cos(alpha_wt_0)
            dw_s = db2_0/np.cos(alpha_wt_0)
            aw0 = 1/2*(dw_s + dw_p0)
            dNf1_0 = np.sqrt((2*aw0*np.sin(alpha_wt_0) - np.sqrt(da2_0**2 - db2_0**2))**2 + db1_0**2) #planeet
            dNf2_0 = np.sqrt((2*aw0*np.sin(alpha_wt_0) - np.sqrt(da1_0**2 - db1_0**2))**2 + db2_0**2) #sun
            a_0 = ((Possible_solutions_pregearing[i][2]+Possible_solutions_pregearing[i][3])/(2*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) + y_max_0)*m_t_0*10**-3*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)
            ksi_Nfw1_0 = np.tan(alpha_wt_0) - np.tan(np.arccos(db1_0/dNf1_0)) #planeet
            ksi_Nfw2_0 = np.tan(alpha_wt_0) - np.tan(np.arccos(db2_0/dNf2_0)) #ring
            ksiNaw1_0 = ksi_Nfw2_0*(Possible_solutions_pregearing[i][2]/Possible_solutions_pregearing[i][3])
            Eps_alpha_0 = (ksi_Nfw1_0 + ksiNaw1_0)/tau_0
            if Eps_alpha_0 < 1 or np.isnan(Eps_alpha_0):
                x_s +=0.005
                if x_s > max_sum_of_prof_shift_0:
                    if len(possibilities_0)==0:
                        bending_stresses_after_PS_0.append([0,0,0,0,0,0,0,0,0,0,0])
                        exit_inner_loops = True
                    else:
                        Best_config_0 = max(possibilities_0, key=lambda x:(x[0]))
                        bending_stresses_after_PS_0.append(Best_config_0)
                        possibilities_0.clear()
                        exit_inner_loops = True
                        exit_inner_loops_second_time = False
                    break
                continue
            Y_F_s,rho_F_s, sFn_s, hFe_s  = calculate_Y_F_ISO('ext',Possible_solutions_pregearing[i][2], Possible_solutions_pregearing[i][9], XE_s_ext, m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_0*10**-3*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180), Eps_alpha_0, x_s)
            Y_F_p0,rho_F_p0, sFn_p0, hFe_p0  = calculate_Y_F_ISO('ext',Possible_solutions_pregearing[i][3], Possible_solutions_pregearing[i][9], XE_P0_ext, m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*10**-3, pressure_angle,0.38*m_t_0*10**-3*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180), Eps_alpha_0, x_P0)
            sigma_s = Possible_solutions_pregearing[i][17]/NOP_pre/(Possible_solutions_pregearing[i][7]*10**-3*m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*10**-3)*Y_F_s*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_pregearing[i][9])*calculate_Y_B_ISO('ext', m_t_0*10**-3*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180), 3.5*m_t_0*10**-3*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_s, rho_F_s, hFe_s)*K
            sigma_P0 = Possible_solutions_pregearing[i][18]/(Possible_solutions_pregearing[i][7]*10**-3*m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*10**-3)*Y_F_p0*calculate_Y_S_ISO(sFn_p0, rho_F_p0, hFe_p0)*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_pregearing[i][9])*calculate_Y_B_ISO('ext', m_t_0*10**-3*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180), 3.5*m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)*10**-3)*Y_DT*K
            sigma_FP_s = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP_pre)*Y_delta_rel_T_ISO(sFn_s/(2*rho_F_s), Material_properties[material][5], calculate_Y_S_ISO(sFn_s, rho_F_s, hFe_s), 1)*Y_R_rel_T_iso(3e6*NOP_pre, material)*Y_X/SF_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            sigma_FP_P0 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p0/(2*rho_F_p0), Material_properties[material][5], calculate_Y_S_ISO(sFn_p0, rho_F_p0, hFe_p0), 1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            sigma_FG_s = sigma_FP_s*SF_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            sigma_FG_P0 = sigma_FP_P0*SF_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            sigma_H_s = calculate_ZH_ISO(Possible_solutions_pregearing[i][9], alpha_wt_0, pressure_angle, 'ext')*Z_E*calculate_Z_eps_ISO(Possible_solutions_pregearing[i][9], Eps_alpha_0)*calculate_Z_beta_ISO(Possible_solutions_pregearing[i][9])*np.sqrt((Possible_solutions_pregearing[i][17]/NOP_pre)/(m_t_0*10**(-3)*Possible_solutions_pregearing[i][2]*Possible_solutions_pregearing[i][7]*10**(-3))*(Possible_solutions_pregearing[i][3]/Possible_solutions_pregearing[i][2]+1)/(Possible_solutions_pregearing[i][3]/Possible_solutions_pregearing[i][2]))*np.sqrt(K)*calculate_Z_B_ISO(alpha_wt_0*180/np.pi,Possible_solutions_pregearing[i][9], Possible_solutions_pregearing[i][2], Possible_solutions_pregearing[i][3], m_t_0*10**(-3), x_s, x_P0, Eps_alpha_0, 'ext')
            sigma_H_P0 = calculate_ZH_ISO(Possible_solutions_pregearing[i][9], alpha_wt_0, pressure_angle, 'ext')*Z_E*calculate_Z_eps_ISO(Possible_solutions_pregearing[i][9], Eps_alpha_0)*calculate_Z_beta_ISO(Possible_solutions_pregearing[i][9])*np.sqrt((Possible_solutions_pregearing[i][17]/NOP_pre)/(m_t_0*10**(-3)*Possible_solutions_pregearing[i][2]*Possible_solutions_pregearing[i][7]*10**(-3))*(Possible_solutions_pregearing[i][3]/Possible_solutions_pregearing[i][2]+1)/(Possible_solutions_pregearing[i][3]/Possible_solutions_pregearing[i][2]))*np.sqrt(K)*calculate_Z_D_ISO(alpha_wt_0, Possible_solutions_pregearing[i][9], Possible_solutions_pregearing[i][2], Possible_solutions_pregearing[i][3],m_t_0*10**(-3), x_P0, Eps_alpha_0)  
            
            sigma_HP_s = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP_pre)/SH_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_pregearing[i][3],Possible_solutions_pregearing[i][2], material, m_t_0*10**(-3), m_t_0*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)), alpha_wt_0)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[i][5])
            sigma_HP_P0 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))*calculate_Z_R_ISO(Possible_solutions_pregearing[i][3],Possible_solutions_pregearing[i][2], material, m_t_0*10**(-3), m_t_0*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)), alpha_wt_0)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[i][5])
            sigma_HG_s = sigma_HP_s*SH_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            sigma_HG_p0 = sigma_HP_P0*SH_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180))
            if sigma_FG_s/(sigma_s*10**(-6)) > SF_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) and  sigma_FG_P0/(sigma_P0*10**(-6)) > SF_min(m_t_0*np.cos(Possible_solutions_pregearing[i][9]*np.pi/180)) and sigma_FG_s/(sigma_s*10**(-6)) > Possible_solutions_pregearing[i][29]:
                if sigma_FG_s/(sigma_s*10**(-6)) -  sigma_FG_P0/(sigma_P0*10**(-6)) <= 0.05 and sigma_FG_s/(sigma_s*10**(-6)) -  sigma_FG_P0/(sigma_P0*10**(-6)) > 0: 
                    possibilities_0.append([(sigma_FG_s/(sigma_s*10**(-6))), sigma_FG_P0/(sigma_P0*10**(-6)), sigma_P0*10**(-6),sigma_s*10**(-6), x_P0, x_s, a_0, sigma_H_s*10**(-6), sigma_H_P0*10**(-6), sigma_HG_s/(sigma_H_s*10**(-6)), sigma_HG_p0/(sigma_H_P0*10**(-6))])
                    x_s +=0.005
                    if x_s > max_sum_of_prof_shift_0:
                        if len(possibilities_0)==0:
                            bending_stresses_after_PS_0.append([0,0,0,0,0,0,0,0,0,0,0])
                            exit_inner_loops = True
                        else:
                            Best_config_0 = max(possibilities_0, key=lambda x:(x[0]))
                            bending_stresses_after_PS_0.append(Best_config_0)
                            possibilities_0.clear()
                            exit_inner_loops = True
                            exit_inner_loops_second_time = False
                        break
                    continue
                else:
                    x_s +=0.005
                    if x_s > max_sum_of_prof_shift_0:
                        if len(possibilities_0)==0:
                            bending_stresses_after_PS_0.append([0,0,0,0,0,0,0,0,0,0,0])
                            exit_inner_loops = True
                        else:
                            Best_config_0 = max(possibilities_0, key=lambda x:(x[0]))
                            bending_stresses_after_PS_0.append(Best_config_0)
                            possibilities_0.clear()
                            exit_inner_loops = True
                            exit_inner_loops_second_time = False
                        break
                    continue
            else:
                x_s +=0.005
                if x_s > max_sum_of_prof_shift_0:
                    if len(possibilities_0)==0:
                        bending_stresses_after_PS_0.append([0,0,0,0,0,0,0,0,0,0,0])
                        exit_inner_loops = True
                    else:
                        Best_config_0 = max(possibilities_0, key=lambda x:(x[0]))
                        bending_stresses_after_PS_0.append(Best_config_0)
                        possibilities_0.clear()
                        exit_inner_loops = True
                        exit_inner_loops_second_time = False
                    break
                continue
        if exit_inner_loops:
            if bending_stresses_after_PS_0[i] == [0, 0, 0, 0, 0, 0 , 0,0, 0, 0, 0]:
                a_0 = Possible_solutions_pregearing[i][13]
            else:
                a_0 = Best_config_0[6]
            #zoeken naar a original in stage 0 
            y__1 = abs(a_0-Possible_solutions_pregearing[i][14])/(m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
            if a_0 < Possible_solutions_pregearing[i][14]:
                y__1 = -y__1
            a_1 = ((Possible_solutions_pregearing[i][0]-Possible_solutions_pregearing[i][1])/(2*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) + y__1)*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)
            cos_alphawt__1 = np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)))/((2*y__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)/(Possible_solutions_pregearing[i][0]-Possible_solutions_pregearing[i][1])) + 1)
            alpha_wt__1 = np.arccos(cos_alphawt__1)
            inv_alpha_wt_min__1 = np.tan(alpha_wt__1) - alpha_wt__1
            sum_of_prof_shift__1 = (inv_alpha_wt_min__1-np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) + np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)))/(2*np.tan(pressure_angle*np.pi/180))*(Possible_solutions_pregearing[i][0] - Possible_solutions_pregearing[i][1]) 
            x_P_1 = -sum_of_prof_shift__1
            step = 0.005
            print(sum_of_prof_shift__1)
            if sum_of_prof_shift__1 < 0:
                while x_P_1 >= sum_of_prof_shift__1:
                    x_R_1 = -sum_of_prof_shift__1 - x_P_1
                    if x_P_1 < 1-(Possible_solutions_pregearing[i][1]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)):
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
                    XE__1_ext = x_P_1 + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    XE__1_int = x_R_1 + Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    da2__1 = m_t__1*10**(-3)*Possible_solutions_pregearing[i][0] - 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1+x_R_1)
                    da1__1 = m_t__1*10**(-3)*Possible_solutions_pregearing[i][1] + 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1+x_P_1)
                    dw_p_1 = db1__1/np.cos(alpha_wt__1) #planet
                    dw_ring_1 = db2__1/np.cos(alpha_wt__1) #ring
                    aw__1 = 1/2*(dw_ring_1 - dw_p_1)
                    dNf1__1 = np.sqrt((2*aw__1*np.sin(alpha_wt__1) - np.sqrt(da2__1**2 - db2__1**2))**2 + db1__1**2)
                    dNf2__1 = np.sqrt((2*aw__1*np.sin(alpha_wt__1) + np.sqrt(da1__1**2 - db1__1**2))**2 + db2__1**2)
                    ksi_Nfw1 = np.tan(alpha_wt__1) - np.tan(np.arccos(db1__1/dNf1__1))
                    ksi_Nfw2 = np.tan(alpha_wt__1) - np.tan(np.arccos(db2__1/dNf2__1))
                    ksiNaw1 = ksi_Nfw2*(-Possible_solutions_pregearing[i][0]/Possible_solutions_pregearing[i][1])
                    Eps_alpha__1 = (ksi_Nfw1 + ksiNaw1)/tau__1
                    alpha_ta1 = np.arccos(db1__1/da1__1) #planeet
                    alpha_ta2 = np.arccos(db2__1/da2__1) #ring
                    beta_a_planeet = np.arctan((np.tan(Possible_solutions_pregearing[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))))/np.cos(alpha_ta1))
                    beta_a_ring = np.arctan((np.tan(Possible_solutions_pregearing[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))))/np.cos(alpha_ta2))
                    d_root_ring =  m_t_1*10**(-3)*Possible_solutions_pregearing[i][0] + 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1.25 - x_R_1)
                    d_root_planet =  m_t_1*10**(-3)*Possible_solutions_pregearing[i][1] - 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1.25 - x_P_1)
                    if Eps_alpha__1 <1 or np.isnan(Eps_alpha__1):
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
                    if a_1 + da1__1/2 > d_root_ring/2 - 0.2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180) or da2__1/2 < d_root_planet/2 + a_1 +0.2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180) or db2__1 >= da2__1:
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
                    if da1__1*(np.pi/(2*Possible_solutions_pregearing[i][1]) + 2*x_P_1*np.tan(pressure_angle*np.pi/180)/Possible_solutions_pregearing[i][1] + np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) - np.tan(alpha_ta1)+ alpha_ta1)*np.cos(beta_a_planeet) < 0.4*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180):
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
                    if da2__1*(np.pi/(2*Possible_solutions_pregearing[i][0]) + 2*x_R_1*np.tan(pressure_angle*np.pi/180)/Possible_solutions_pregearing[i][0] - (np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) - np.tan(alpha_ta2)+ alpha_ta2))*np.cos(beta_a_ring) < 0.4*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180):
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
                    Y_F_R_1,rho_F_R_1, sFn_R_1, hFe_R_1  = calculate_Y_F_ISO('int',Possible_solutions_pregearing[i][0], Possible_solutions_pregearing[i][10], XE__1_int, m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), Eps_alpha__1, x_R_1)
                    Y_F_p_1,rho_F_p_1, sFn_p_1, hFe_p_1  = calculate_Y_F_ISO('ext',Possible_solutions_pregearing[i][1], Possible_solutions_pregearing[i][10], XE__1_ext, m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), Eps_alpha__1, x_P_1)
                    sigma_R_1 = Possible_solutions_pregearing[i][15]/NOP_pre/(Possible_solutions_pregearing[i][8]*10**-3*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3)*Y_F_R_1*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_pregearing[i][10])*calculate_Y_B_ISO('int', m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), 3.5*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R_1, rho_F_R_1, hFe_R_1)*K
                    sigma_p_1 = Possible_solutions_pregearing[i][16]/(Possible_solutions_pregearing[i][8]*10**-3*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3)*Y_F_p_1*calculate_Y_S_ISO(sFn_p_1, rho_F_p_1, hFe_p_1)*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_pregearing[i][10])*calculate_Y_B_ISO('ext', m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), 3.5*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3)*Y_DT*K
                    sigma_FP_R_1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP_pre)*Y_delta_rel_T_ISO(sFn_R_1/(2*rho_F_R_1), Material_properties[material][5], calculate_Y_S_ISO(sFn_R_1, rho_F_R_1, hFe_R_1), 1)*Y_R_rel_T_iso(3e6*NOP_pre, material)*Y_X/SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_FP_P_1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p_1/(2*rho_F_p_1), Material_properties[material][5], calculate_Y_S_ISO(sFn_p_1, rho_F_p_1, hFe_p_1), 1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_FG_R_1 = sigma_FP_R_1*SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_FG_p_1 = sigma_FP_P_1*SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))

                    sigma_H_R_1 = calculate_ZH_ISO(Possible_solutions_pregearing[i][10], alpha_wt__1, pressure_angle, 'int')*Z_E*calculate_Z_eps_ISO(Possible_solutions_pregearing[i][10], Eps_alpha__1)*calculate_Z_beta_ISO(Possible_solutions_pregearing[i][10])*np.sqrt((Possible_solutions_pregearing[i][15]/NOP_pre)/(m_t__1*10**(-3)*Possible_solutions_pregearing[i][1]*Possible_solutions_pregearing[i][8]*10**-3)*(-Possible_solutions_pregearing[i][0]/Possible_solutions_pregearing[i][1]+1)/(-Possible_solutions_pregearing[i][0]/Possible_solutions_pregearing[i][1]))*np.sqrt(K) #Z_D is = 1 for internal gears
                    sigma_H_P_1 = sigma_H_R_1*calculate_Z_B_ISO(alpha_wt__1*180/np.pi, Possible_solutions_pregearing[i][10], Possible_solutions_pregearing[i][1], -Possible_solutions_pregearing[i][0],m_t__1*10**(-3), x_P_1, x_R_1, Eps_alpha__1, 'int')
                    sigma_HP_R_1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP_pre)/SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))*calculate_Z_R_ISO( Possible_solutions_pregearing[i][1],-Possible_solutions_pregearing[i][0], material, m_t__1*10**(-3), m_t__1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)), alpha_wt__1)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[i][6])
                    sigma_HP_P_1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))*calculate_Z_R_ISO( Possible_solutions_pregearing[i][1],-Possible_solutions_pregearing[i][0], material, m_t__1*10**(-3), m_t__1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)), alpha_wt__1)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[i][6])
                    sigma_HG_R_1 = sigma_HP_R_1*SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_HG_P_1 = sigma_HP_P_1*SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    if sigma_FG_R_1/(sigma_R_1*10**(-6)) > SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) and  sigma_FG_p_1/(sigma_p_1*10**(-6)) > SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)):
                        possibilities__1.append([(sigma_FG_R_1/(sigma_R_1*10**(-6))), sigma_FG_p_1/(sigma_p_1*10**(-6)), sigma_p_1*10**(-6), sigma_R_1*10**(-6), x_P_1, x_R_1, a_1, sigma_H_R_1*10**(-6), sigma_H_P_1*10**(-6), sigma_HP_R_1, sigma_HP_P_1, sigma_HG_R_1/(sigma_H_R_1*10**(-6)), sigma_HG_P_1/(sigma_H_P_1*10**(-6))])
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
                    else:
                        x_P_1 -= step
                        if x_P_1 < sum_of_prof_shift__1:
                            if len(possibilities__1)==0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = max(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                            break
                        continue
            else:
                step = -0.005
                x_P_1 = -sum_of_prof_shift__1
                if x_P_1 <1-(Possible_solutions_pregearing[i][1]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)):
                    x_P_1 = 1-(Possible_solutions_pregearing[i][1]*np.sin(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)))**2)/(2*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    # bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                    # exit_inner_loops_second_time = True
                while x_P_1 <= sum_of_prof_shift__1:
                    x_R_1 = -sum_of_prof_shift__1 - x_P_1
                    XE__1_ext = x_P_1 + Tooth_thickness_allowances[allowances][1]/(2*np.tan(pressure_angle*np.pi/180)*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    XE__1_int = x_R_1 + Tooth_thickness_allowances[allowances][3]/(2*np.tan(pressure_angle*np.pi/180)*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    da2__1 = m_t__1*10**(-3)*Possible_solutions_pregearing[i][0] - 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1+x_R_1)
                    da1__1 = m_t__1*10**(-3)*Possible_solutions_pregearing[i][1] + 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1+x_P_1)
                    dw_p_1 = db1__1/np.cos(alpha_wt__1)
                    dw_ring_1 = db2__1/np.cos(alpha_wt__1)
                    aw__1 = 1/2*(dw_ring_1 - dw_p_1)
                    dNf1__1 = np.sqrt((2*aw__1*np.sin(alpha_wt__1) - np.sqrt(da2__1**2 - db2__1**2))**2 + db1__1**2)
                    dNf2__1 = np.sqrt((2*aw__1*np.sin(alpha_wt__1) + np.sqrt(da1__1**2 - db1__1**2))**2 + db2__1**2)
                    ksi_Nfw1 = np.tan(alpha_wt__1) - np.tan(np.arccos(db1__1/dNf1__1))
                    ksi_Nfw2 = np.tan(alpha_wt__1) - np.tan(np.arccos(db2__1/dNf2__1))
                    ksiNaw1 = ksi_Nfw2*(-Possible_solutions_pregearing[i][0]/Possible_solutions_pregearing[i][1])
                    Eps_alpha__1 = (ksi_Nfw1 + ksiNaw1)/tau__1
                    if Eps_alpha__1 <1 or np.isnan(Eps_alpha__1):
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
                    alpha_ta1 = np.arccos(db1__1/da1__1) #planeet
                    alpha_ta2 = np.arccos(db2__1/da2__1) #ring
                    beta_a_planeet = np.arctan((np.tan(Possible_solutions_pregearing[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))))/np.cos(alpha_ta1))
                    beta_a_ring = np.arctan((np.tan(Possible_solutions_pregearing[i][10]*np.pi/180) * np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))))/np.cos(alpha_ta2))
                    d_root_ring =  m_t__1*10**(-3)*Possible_solutions_pregearing[i][0] + 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1.25 - x_R_1)
                    d_root_planet =  m_t__1*10**(-3)*Possible_solutions_pregearing[i][1] - 2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*(1.25 - x_P_1)
                    if Eps_alpha__1 <1 or np.isnan(Eps_alpha__1):
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
                    if a_1 + da1__1/2 > d_root_ring/2 - 0.2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180) or da2__1/2 < d_root_planet/2 + a_1 +0.2*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180) or db2__1 >= da2__1:
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
                    if da1__1*(np.pi/(2*Possible_solutions_pregearing[i][1]) + 2*x_P_1*np.tan(pressure_angle*np.pi/180)/Possible_solutions_pregearing[i][1] + np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) - np.tan(alpha_ta1)+ alpha_ta1)*np.cos(beta_a_planeet) < 0.4*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180):
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
                    if da2__1*(np.pi/(2*Possible_solutions_pregearing[i][0]) + 2*x_R_1*np.tan(pressure_angle*np.pi/180)/Possible_solutions_pregearing[i][0] - (np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))) - np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) - np.tan(alpha_ta2)+ alpha_ta2))*np.cos(beta_a_ring) < 0.4*m_t__1*10**(-3)*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180):
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
                    Y_F_R_1,rho_F_R_1, sFn_R_1, hFe_R_1  = calculate_Y_F_ISO('int',Possible_solutions_pregearing[i][0], Possible_solutions_pregearing[i][10], XE__1_int, m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), Eps_alpha__1, x_R_1)
                    Y_F_p_1,rho_F_p_1, sFn_p_1, hFe_p_1  = calculate_Y_F_ISO('ext',Possible_solutions_pregearing[i][1], Possible_solutions_pregearing[i][10], XE__1_ext, m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3, pressure_angle,0.38*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), Eps_alpha__1, x_P_1)
                    sigma_R_1 = Possible_solutions_pregearing[i][15]/NOP_pre/(Possible_solutions_pregearing[i][8]*10**-3*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3)*Y_F_R_1*calculate_Y_beta_ISO(overlap_ratio, Possible_solutions_pregearing[i][10])*calculate_Y_B_ISO('int', m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), 3.5*m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))*Y_DT*calculate_Y_S_ISO(sFn_R_1, rho_F_R_1, hFe_R_1)*K
                    sigma_p_1 = Possible_solutions_pregearing[i][16]/(Possible_solutions_pregearing[i][8]*10**-3*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3)*Y_F_p_1*calculate_Y_S_ISO(sFn_p_1, rho_F_p_1, hFe_p_1)*calculate_Y_beta_ISO( overlap_ratio, Possible_solutions_pregearing[i][10])*calculate_Y_B_ISO('ext', m_t__1*10**-3*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180), 3.5*m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)*10**-3)*Y_DT*K
                    sigma_FP_R_1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6*NOP_pre)*Y_delta_rel_T_ISO(sFn_R_1/(2*rho_F_R_1), Material_properties[material][5], calculate_Y_S_ISO(sFn_R_1, rho_F_R_1, hFe_R_1), 1)*Y_R_rel_T_iso(3e6*NOP_pre, material)*Y_X/SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_FP_P_1 = Material_properties[material][4]*Y_ST*calculate_Y_NT_ISO(3*10**6)*Y_delta_rel_T_ISO(sFn_p_1/(2*rho_F_p_1), Material_properties[material][5], calculate_Y_S_ISO(sFn_p_1, rho_F_p_1, hFe_p_1), 1)*Y_R_rel_T_iso(3e6, material)*Y_X/SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_FG_R_1 = sigma_FP_R_1*SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_FG_p_1 = sigma_FP_P_1*SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))

                    sigma_H_R_1 = calculate_ZH_ISO(Possible_solutions_pregearing[i][10], alpha_wt__1, pressure_angle, 'int')*Z_E*calculate_Z_eps_ISO(Possible_solutions_pregearing[i][10], Eps_alpha__1)*calculate_Z_beta_ISO(Possible_solutions_pregearing[i][10])*np.sqrt((Possible_solutions_pregearing[i][15]/NOP_pre)/(m_t__1*10**(-3)*Possible_solutions_pregearing[i][1]*Possible_solutions_pregearing[i][8]*10**-3)*(-Possible_solutions_pregearing[i][0]/Possible_solutions_pregearing[i][1]+1)/(-Possible_solutions_pregearing[i][0]/Possible_solutions_pregearing[i][1]))*np.sqrt(K) #Z_D is = 1 for internal gears
                    sigma_H_P_1 = sigma_H_R_1*calculate_Z_B_ISO(alpha_wt__1*180/np.pi, Possible_solutions_pregearing[i][10], Possible_solutions_pregearing[i][1], -Possible_solutions_pregearing[i][0],m_t__1*10**(-3), x_P_1, x_R_1, Eps_alpha__1, 'int')
                    sigma_HP_R_1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9*NOP_pre)/SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))*calculate_Z_R_ISO( Possible_solutions_pregearing[i][1],-Possible_solutions_pregearing[i][0], material, m_t__1*10**(-3), m_t__1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)), alpha_wt__1)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[i][6])
                    sigma_HP_P_1 = Material_properties[material][7]*calculate_Z_NT_ISO(10**9)/SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))*calculate_Z_R_ISO( Possible_solutions_pregearing[i][1],-Possible_solutions_pregearing[i][0], material, m_t__1*10**(-3), m_t__1*10**(-3), np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)), alpha_wt__1)*Z_W*Z_X*calculate_ZL_ISO(material)*calculate_Z_V_ISO(speeds_pre[i][6])
                    sigma_HG_R_1 = sigma_HP_R_1*SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    sigma_HG_P_1 = sigma_HP_P_1*SH_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180))
                    if sigma_FG_R_1/(sigma_R_1*10**(-6)) > SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)) and  sigma_FG_p_1/(sigma_p_1*10**(-6)) > SF_min(m_t__1*np.cos(Possible_solutions_pregearing[i][10]*np.pi/180)):
                        possibilities__1.append([(sigma_FG_R_1/(sigma_R_1*10**(-6))), sigma_FG_p_1/(sigma_p_1*10**(-6)), sigma_p_1*10**(-6), sigma_R_1*10**(-6), x_P_1, x_R_1, a_1, sigma_H_R_1*10**(-6), sigma_H_P_1*10**(-6), sigma_HP_R_1, sigma_HP_P_1, sigma_HG_R_1/(sigma_H_R_1*10**(-6)), sigma_HG_P_1/(sigma_H_P_1*10**(-6))])
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
                    else:
                        x_P_1 -= step
                        if x_P_1 > sum_of_prof_shift__1:
                            if len(possibilities__1) == 0:
                                bending_stresses_after_PS_0[i].extend([0,0,0,0,0,0,0,0,0,0,0,0,0])
                                exit_inner_loops_second_time = True
                            else:
                                Best_config__1 = min(possibilities__1, key=lambda x:(abs(x[0]-x[1])))
                                bending_stresses_after_PS_0[i].extend(Best_config__1)
                                possibilities__1.clear()
                                exit_inner_loops_second_time = True
                            break
                        continue
        if exit_inner_loops_second_time:
            continue
    return bending_stresses_after_PS_0

def calculate_pre_gearing():
    #calculating the pregearing
    global sorted_final_solutions_pregearing, wolfrom_configuration_chosen
    values = [entry.get() for entry in entry_boxes]
    if module0ok:
        m0 = float(values[14]) #is m_t in het geval van helical gears0
    else:
        if float(values[14]) <= m_min_bending_stage_0*10**3:
            m0 = round(m_min_bending_stage_0*10**3,4)
        else:
            m0 = float(values[7]) #is m_t in het geval van helical gears

    if module_1ok:
        m_1 = float(values[15])#is m_t in het geval van helical gears
    else:
        if float(values[15]) <= m_min_bending_stage__1*10**3:
            m_1 = round(m_min_bending_stage__1*10**3,4)
        else:
            m_1 = float(values[15]) 
    Desired_GR_pregearing = int(values[11])
    Deviation = int(values[12])
    NOP_pre = int(values[13])
    sin_pre = np.sin(np.pi/NOP_pre)
    m_1 = float(values[15]) #module of 1' stage
    Final_solutions_pregearing = []
    pregearing_setpoint = float(values[16])/100
    m0 = float(values[14])
    WR2 = float(values[2])
    max_diam = int(values[0])
    axial_length = float(values[10])
    material = material_var.get()
    allowances = allowance_var.get()
    Type_of_gear = gear_var.get()
    m_t__1 = m_1
    m_t_0 = m0
    hole_decision = hole_var.get()
    if hole_decision =='Yes':
        hole_diameter = int(values[19])
        wolfrom_configuration_chosen = sorted_possible_solutions_wolfrom[int(values[17])-1]
    else:
        hole_diameter = 0
        wolfrom_configuration_chosen = sorted_possible_solutions_wolfrom[int(values[17])-1]

    eff_w = wolfrom_configuration_chosen[5]
    GR_w = wolfrom_configuration_chosen[4]
    ZR2 = wolfrom_configuration_chosen[0]
    ZR1 = wolfrom_configuration_chosen[1]
    Tr2 = float(values[1])

    for Zr_1 in range(int(max_diam/m_t__1) ,22,-1):
        diameters = [ZR2*m_t_2, ZR1*m_t_1, Zr_1*m_t__1]
        max_diam_ = max(diameters)
        max_index = diameters.index(max_diam_)
        if max_index == 0:
            m_max = m_t_2
        elif max_index == 1:
            m_max = m_t_1
        else:
            m_max = m_t__1
        for Zp_1 in range(int(Zr_1/2 - hole_diameter/(2*m_t__1)),11,-1):
            for Zp0 in range(int(((max_diam_/2 + 3.5*m_max + 1.25*m_max*np.cos(0) - m_t_0)/m_t_0)-11/2), 11, -1):
                for Zs in range(int((max_diam_ + 2*3.5*m_max + 2*1.25*m_max*np.cos(0) - 2*m_t_0)/m_t_0)-22, 11, -1):
                    S_1 = Zr_1/Zp_1
                    S0 = Zp0/Zs
                    if (Desired_GR_pregearing - Desired_GR_pregearing*Deviation/100) <= (1+S0*S_1) <= (Desired_GR_pregearing + Desired_GR_pregearing*Deviation/100):
                        if Type_of_gear == 'Spur':
                            overlap_ratio = 0
                            Helix_angle_input = 0           
                            Helix_angle_output = 0
                            beta_b__1 = np.arctan(np.tan(Helix_angle_output*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output*np.pi/180))))
                            beta_b_0 = np.arctan(np.tan(Helix_angle_input*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input*np.pi/180))))
                            p_b_0_t = np.pi*m_t_0*10**(-3)*np.cos(pressure_angle*np.pi/180)
                            p_b_1__t = np.pi*m_t__1*10**(-3)*np.cos(pressure_angle*np.pi/180)
                            da2__1 = m_t__1*10**(-3)*Zr_1 - 2*m_t__1*10**(-3) #for the internal one
                            db2__1 = m_t__1*10**(-3)*Zr_1*np.cos(pressure_angle*math.pi/180) #for the internal one
                            da1__1 = m_t__1*10**(-3)*Zp_1 + 2*m_t__1*10**(-3) #for the external one
                            db1__1 = m_t__1*10**(-3)*Zp_1*np.cos(pressure_angle*math.pi/180) #for the external one
                            a__1 = (m_t__1*10**(-3)*Zr_1 - m_t__1*10**(-3)*Zp_1)/2 
                            # #transverse contact ratio input mesh
                            da2_0 = m_t_0*10**(-3)*Zp0 + 2*m_t_0*10**(-3) # 1 is sun, 2 is planet
                            db2_0= m_t_0*10**(-3)*Zp0*np.cos(pressure_angle*math.pi/180)
                            da1_0 = m_t_0*10**(-3)*Zs + 2*m_t_0*10**(-3)
                            db1_0 = m_t_0*10**(-3)*Zs*np.cos(pressure_angle*math.pi/180) 
                            a_0 = (m_t_0*10**(-3)*Zp0 + m_t_0*10**(-3)*Zs)/2 
                            Eps_alpha_0 = 1/p_b_0_t*(np.sqrt(da1_0**2 - db1_0**2)/2 + np.sqrt(da2_0**2 - db2_0**2)/2 - a_0*np.sin(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input*np.pi/180))))
                            Eps_alpha__1 = 1/p_b_1__t*(np.sqrt(da1__1**2-db1__1**2)/2 - np.sqrt(da2__1**2-db2__1**2)/2 + a__1*np.sin(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output*np.pi/180))))
                            Epsilon_approach_0 = Zs/(2*np.pi) * (np.tan(np.arccos(db1_0/da1_0)) - np.tan(pressure_angle*np.pi/180))
                            Epsilon_recess_0 = Zp0/(2*np.pi) * (np.tan(np.arccos(db2_0/da2_0)) - np.tan(pressure_angle*np.pi/180))
                            Epsilon_approach__1 = -Zr_1/(2*np.pi) * (np.tan(np.arccos(db2__1/da2__1)) - np.tan(pressure_angle*np.pi/180))
                            Epsilon_recess__1 = Zp_1/(2*np.pi) * (np.tan(np.arccos(db1__1/da1__1)) - np.tan(pressure_angle*np.pi/180))
                            f_R_1_P_1 = fric_coeff*np.pi*(1/Zp_1 - 1/Zr_1)*(1-Eps_alpha__1 + Epsilon_approach__1**2 + Epsilon_recess__1**2)/np.cos(beta_b__1)
                            f_S_P0 = fric_coeff*np.pi*(1/Zp0 + 1/Zs)*(1-Eps_alpha_0 + Epsilon_approach_0**2 + Epsilon_recess_0**2)/np.cos(beta_b_0)
                            eff_real = ((1-(1-1/(1+S0*S_1))*f_S_P0 + 1/(1+S0*S_1)*f_R_1_P_1)/(f_R_1_P_1+1))
                        else:
                            if Type_of_gear == "Helical, eps_beta = 1":
                                overlap_ratio = 1
                            if Type_of_gear == "Helical, eps_beta = 2":
                                overlap_ratio = 2
                            error_1 = 50
                            error_2 = 50
                            Helix_angle_input = 0           
                            Helix_angle_output = 0
                            a__1 = (m_t__1*10**(-3)*Zr_1 - m_t__1*10**(-3)*Zp_1)/2 
                            a_0 = (m_t_0*10**(-3)*Zp0 + m_t_0*10**(-3)*Zs)/2 
                            beta_b__1_start = np.arctan(np.tan(Helix_angle_output*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output*np.pi/180))))
                            beta_b_0_start = np.arctan(np.tan(Helix_angle_input*np.pi/180)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input*np.pi/180))))
                            while error_1 and error_2 >10**(-8):
                                da2__1 = m_t__1*10**(-3)*Zr_1 - 2*m_t__1*10**(-3)*np.cos(Helix_angle_output) #for the internal one
                                da1__1 = m_t__1*10**(-3)*Zp_1 + 2*m_t__1*10**(-3)*np.cos(Helix_angle_output) #for the external one
                                da2_0 = m_t_0*10**(-3)*Zp0 + 2*m_t_0*10**(-3)*np.cos(Helix_angle_input) # 1 is sun, 2 is planet
                                da1_0 = m_t_0*10**(-3)*Zs + 2*m_t_0*10**(-3)*np.cos(Helix_angle_input)
                                db2__1 = m_t__1*10**(-3)*Zr_1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))) #for the internal one
                                db1__1 = m_t__1*10**(-3)*Zp_1*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))) #for the external one
                                db2_0 = m_t_0*10**(-3)*Zp0*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input)))
                                db1_0 = m_t_0*10**(-3)*Zs*np.cos(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input)))
                                Epsilon_approach_0 = Zs/(2*np.pi)*(np.tan(np.arccos(db1_0/da1_0)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))))
                                Epsilon_recess_0 = Zp0/(2*np.pi)*(np.tan(np.arccos(db2_0/da2_0)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_input))))
                                Epsilon_approach__1 = -Zr_1/(2*np.pi) * (np.tan(np.arccos(db2__1/da2__1)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))))
                                Epsilon_recess__1 = Zp_1/(2*np.pi) * (np.tan(np.arccos(db1__1/da1__1)) - np.tan(np.arctan(np.tan(pressure_angle*math.pi/180)/np.cos(Helix_angle_output))))
                                p_b_0 = np.pi*m_t_0*10**(-3)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input)))
                                p_b_1_ = np.pi*m_t__1*10**(-3)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output)))
                                Eps_alpha_0 = 1/p_b_0*(np.sqrt(da1_0**2 - db1_0**2)/2 + np.sqrt(da2_0**2 - db2_0**2)/2 - a_0*np.sin(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input))))
                                Eps_alpha__1 = 1/p_b_1_*(np.sqrt(da1__1**2-db1__1**2)/2 - np.sqrt(da2__1**2-db2__1**2)/2 + a__1*np.sin(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output))))
                                f_R_1_P_1 = (fric_coeff*np.pi*(1/Zp_1 - 1/Zr_1)*(1-Eps_alpha__1 + Epsilon_approach__1**2 + Epsilon_recess__1**2))/np.cos(beta_b__1_start)
                                f_S_P0 = (fric_coeff*np.pi*(1/Zp0 + 1/Zs)*(1-Eps_alpha_0 + Epsilon_approach_0**2 + Epsilon_recess_0**2))/np.cos(beta_b_0_start)
                                eff_real = (1-(1-1/(1+S0*S_1))*f_S_P0 + 1/(1+S0*S_1)*f_R_1_P_1)/(f_R_1_P_1+1)
                                GR_0 = (1+S0*S_1)
                                eff_tot = eff_real*eff_w
                                Fs = (-1/eff_tot * 1/(GR_0*GR_w)*Tr2)/(m_t_0*10**-3*Zs/2)
                                Fp0 = Fs/NOP_pre
                                Fr_1 = ((1/(GR_0*GR_w)*1/eff_tot + 1/GR_w*1/eff_w)*Tr2)/(m_t__1*10**-3*Zr_1/2)
                                Fp_1 = Fr_1/NOP_pre
                                b0, b_1 = calculate_width(Fp0, Fp_1, axial_length, Type_of_gear, m_t_0, m_t__1)
                                Helix_angle_input = np.arctan(overlap_ratio*np.pi*m_t_0/b0)         
                                Helix_angle_output = np.arctan(overlap_ratio*np.pi*m_t__1/b_1)
                                beta_b_0 = np.arctan(np.tan(Helix_angle_input)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_input))))
                                beta_b__1 = np.arctan(np.tan(Helix_angle_output)*np.cos(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(Helix_angle_output))))
                                error_1 = abs(beta_b_0 - beta_b_0_start)
                                error_2 = abs(beta_b__1 - beta_b__1_start)
                                beta_b__1_start = beta_b__1
                                beta_b_0_start = beta_b_0
                        fR1_P1__ = fric_coeff*np.pi*(1/Zp_1 - 1/Zr_1)*0.5/np.cos(beta_b__1)
                        fS_P0_ = fric_coeff*np.pi*(1/Zp0 + 1/Zs)*0.5/np.cos(beta_b_0)
                        eff = ((1-(1-1/(1+S0*S_1))*fS_P0_ + 1/(1+S0*S_1)*fR1_P1__)/(fR1_P1__+1))
                        if eff >= pregearing_setpoint:
                            if Eps_alpha__1 < 1 or Eps_alpha_0 <1:
                                continue
                            if Zs*m_t_0 + 2*Zp0*m_t_0 + 2*m_t_0*np.cos(Helix_angle_input) > max_diam_ + 2*3.5*m_max + 2 * 1.25 * m_max*np.cos(Helix_angle_input):
                                continue
                            if abs((Zr_1*m_t__1 - Zp_1*m_t__1) - (Zp0*m_t_0 + Zs*m_t_0)) > min(m_t__1,m_t_0): #coaxiality (in de veronderstelling dat m1=m2)
                                continue
                            if hole_decision == 'Yes':
                                    if (Zr_1*m_t__1/2-2*(Zp_1*m_t__1/2) - m_t__1*np.cos(Helix_angle_output)) < hole_diameter/2 + margin: #clearance/2 want je moet maar langs 1 kant clearance hebben
                                        continue
                                    if Zs*m_t_0/2 - hole_diameter/2 < 3.5*m_t_0: #neighbouring voor cable hole
                                        continue    
                            if (Zs + Zp0)*sin_pre < Zp0 + 2*np.cos(Helix_angle_input) + clearance/m_t_0: #neighbouring
                                continue
                            if (Zr_1 - Zp_1)*sin_pre < Zp_1 + 2*np.cos(Helix_angle_output) + clearance/m_t__1: #neighbouring
                                continue
                            if (Zp0*Zr_1+Zp_1*Zs)/NOP_pre != int((Zp0*Zr_1+Zp_1*Zs)/NOP_pre): #assembly condition, when identical planets are assumed
                                continue
                            else:
                                Final_solutions_pregearing.append([Zr_1, Zp_1, Zs, Zp0, (1+S0*S_1), eff_real, eff]) 
                                Epsilon_alpha_0.append(Eps_alpha_0)
                                Epsilon_alpha__1.append(Eps_alpha__1) 
                                HA__1.append(Helix_angle_output*180/np.pi)
                                HA_0.append(Helix_angle_input*180/np.pi)
                                A_pre.append([a_0,a__1])
                      
    Forces_pregearing = calculate_tangential_force_pregearing(Final_solutions_pregearing,  values, m_t_0, m_t__1, wolfrom_configuration_chosen, NOP_pre)
    bending_stresses_pregearing, permissible_bending_stresses_pregearing, SF_pregearing = calculate_bending_stress_pregearing(Final_solutions_pregearing, m_t_0, m_t__1, Forces_pregearing, NOP_pre, HA_0, HA__1, axial_length, Type_of_gear, material, allowances)
    speeds_pre = calculate_speeds_pregearing(Final_solutions_pregearing, wolfrom_configuration_chosen, WR2, m_t_0, m_t__1)
    Hertzian_stresses_pregearing, permissible_H_pre, safety_H_pre = calculate_Hertzian_Stress_pregearing(Final_solutions_pregearing, m_t_0, m_t__1, Forces_pregearing, NOP_pre, HA_0, HA__1, pressure_angle, material,0, B0, B_1, speeds_pre)
    
    for i in range(len(Forces_pregearing)):
        Final_solutions_pregearing[i].append(B0[i])
        Final_solutions_pregearing[i].append(B_1[i])
        Final_solutions_pregearing[i].append(HA_0[i])
        Final_solutions_pregearing[i].append(HA__1[i])
        Final_solutions_pregearing[i].append(Epsilon_alpha_0[i])
        Final_solutions_pregearing[i].append(Epsilon_alpha__1[i])
        Final_solutions_pregearing[i].append(A_pre[i][0])
        Final_solutions_pregearing[i].append(A_pre[i][1])
        Final_solutions_pregearing[i].append(Forces_pregearing[i][0])
        Final_solutions_pregearing[i].append(Forces_pregearing[i][1])
        Final_solutions_pregearing[i].append(Forces_pregearing[i][2])
        Final_solutions_pregearing[i].append(Forces_pregearing[i][3])
        Final_solutions_pregearing[i].append(bending_stresses_pregearing[i][0])
        Final_solutions_pregearing[i].append(bending_stresses_pregearing[i][1])
        Final_solutions_pregearing[i].append(bending_stresses_pregearing[i][2])
        Final_solutions_pregearing[i].append(bending_stresses_pregearing[i][3])
        Final_solutions_pregearing[i].append(permissible_bending_stresses_pregearing[i][0])
        Final_solutions_pregearing[i].append(permissible_bending_stresses_pregearing[i][1])
        Final_solutions_pregearing[i].append(permissible_bending_stresses_pregearing[i][2])
        Final_solutions_pregearing[i].append(permissible_bending_stresses_pregearing[i][3])
        Final_solutions_pregearing[i].append(SF_pregearing[i][0])
        Final_solutions_pregearing[i].append(SF_pregearing[i][1])
        Final_solutions_pregearing[i].append(SF_pregearing[i][2])
        Final_solutions_pregearing[i].append(SF_pregearing[i][3])
        Final_solutions_pregearing[i].append(Hertzian_stresses_pregearing[i][0])
        Final_solutions_pregearing[i].append(Hertzian_stresses_pregearing[i][1])
        Final_solutions_pregearing[i].append(Hertzian_stresses_pregearing[i][2])
        Final_solutions_pregearing[i].append(Hertzian_stresses_pregearing[i][3])
        Final_solutions_pregearing[i].append(permissible_H_pre[i][0])
        Final_solutions_pregearing[i].append(permissible_H_pre[i][1])
        Final_solutions_pregearing[i].append(permissible_H_pre[i][2])
        Final_solutions_pregearing[i].append(permissible_H_pre[i][3])
        Final_solutions_pregearing[i].append(safety_H_pre[i][0])
        Final_solutions_pregearing[i].append(safety_H_pre[i][1])
        Final_solutions_pregearing[i].append(safety_H_pre[i][2])
        Final_solutions_pregearing[i].append(safety_H_pre[i][3])

    sorted_final_solutions_pregearing = sorted(Final_solutions_pregearing, key=lambda x: x[6], reverse = True)
    sorted_final_solutions_pregearing = sorted_final_solutions_pregearing[:20]
    if len(Final_solutions_pregearing) < 20:
        sorted_final_solutions_pregearing = sorted_final_solutions_pregearing[:len(Final_solutions_pregearing)]
    After_PS = profile_shift_pregearing(material, sorted_final_solutions_pregearing, m_t_0, m_t__1, NOP_pre, allowances, speeds_pre)

    for i in range(len(sorted_final_solutions_pregearing)):
        for j in range(24):
            sorted_final_solutions_pregearing[i].append(After_PS[i][j])

    for i in range(len(sorted_final_solutions_pregearing)):  #43 is eerste
        if sorted_final_solutions_pregearing[i][60] > max_diam*10**(-3)/2 - sorted_final_solutions_pregearing[i][1]*m_t__1*10**(-3)/2:
            sorted_final_solutions_pregearing[i].append('After profile shift, the solution does not fit anymore, this is not a good solution!')
        else: 
            sorted_final_solutions_pregearing[i].append('After profile shift, the solution still fits!')

        if sorted_final_solutions_pregearing[i][60] == 0:
            sorted_final_solutions_pregearing[i].append('This is not a good solution, the profile shift could not be used to comply with the coaxiality condition or to prevent undercutting.')

    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, f"There are {len(Final_solutions_pregearing)} solutions", "bold")
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "Possible configurations for pregearing: ", "bold")
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "\n", "bold")
    header = "{:<7}{:<5}{:<5}{:<5}{:<5}{:<8}{:<15}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}".format(
        "#", "Zr_1", "Zp_1", "Zs", "Zp0", "GR", "Efficiency","eff_ea=1", "b0", 'b_1', 'HA_0', 'HA__1', 'Eps_a_0', 'Eps_a__1','a_0', 'a__1', "Ft_R_1[N]", "Ft_P_1[N]", "Ft_s[N]", "Ft_p0[N]",
        "BS_R_1[MPa]", "BS_P_1[MPa]", "BS_s[MPa]", "BS_P0[MPa]", "BS_R_1_perm", "BS_P_1_perm", "BS_s_perm", "BS_P0_perm", "SF_R_1", "SF_P_1", "SF_s", "SF_P0","HS_R_1[MPa]", "HS_P_1[MPa]", "HS_s[MPa]", "HS_P0[MPa]", 'HS_R_1_perm',
          "HS_P_1_perm", "HS_s_perm", "HS_P0_perm", "SF_H_R_1", "SF_H_P_1", "SF_H_s", "SF_H_P0",  'SF_x_s', 'SF_x_P0', 'BS_P0', 'BS_s', 'x_P0', 'x_s', 'a_0', 'sigma_H_s', 'sigma_H_P0', 'S_H_x_s', 'S_H_x_P0', 'SF_x_R_1', 'SF_x_P_1', 'BS_P_1', 'BS_R_1', 'x_P_1', 'x_R_1', 'a__1', 'sigma_H_R_1', 'sigma_H_P_1','sigma_HP_R1', 'sigma_HP_P1', 'S_H_x_R_1', 'S_H_x_P_1', 'Warning' )
    output_text.insert(tk.END, header + "\n")

    output_text.tag_configure("red", foreground="red")
                                                                                                                                                                                                                          
    for solution in sorted_final_solutions_pregearing[:]:
        result = "{:<7}{:<5}{:<5}{:<5}{:<5}{:<8}{:<15}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{:<12}{}".format(
            sorted_final_solutions_pregearing.index(solution)+1, solution[0], solution[1], solution[2], solution[3], round(solution[4], 4), round(solution[5], 4),
            round(solution[6], 4), round(solution[7], 4), round(solution[8], 4), round(solution[9], 4),
            round(solution[10], 4), round(solution[11], 4), round(solution[12], 3), round(solution[13], 3),
            round(solution[14], 3), round(solution[15], 3), round(solution[16], 3), round(solution[17], 3), 
            round(solution[18],3), round(solution[19], 3), round(solution[20], 3), round(solution[21], 3), 
            round(solution[22], 3), round(solution[23], 3), round(solution[24],3), round(solution[25],3), 
            round(solution[26],3), round(solution[27],3), round(solution[28],3), round(solution[29],3), 
            round(solution[30],3), round(solution[31],3), round(solution[32],3), round(solution[33],3), round(solution[34],3),
            round(solution[35],3), round(solution[36],3), round(solution[37],3), round(solution[38],3), round(solution[39],3), round(solution[40],3), round(solution[41],3), round(solution[42],3),
            round(solution[43],3), round(solution[44],3), round(solution[45],3), round(solution[46],3), round(solution[47],3), round(solution[48],3), round(solution[49],3), round(solution[50],3),
            round(solution[51],3), round(solution[52],3), round(solution[53],3), round(solution[54],3), round(solution[55],3), round(solution[56],3), round(solution[57],3), round(solution[58],3), 
            round(solution[59],3), round(solution[60],3), round(solution[61],3), round(solution[62],3), round(solution[63],3), round(solution[64],3), round(solution[65],3), round(solution[66],3), solution[67])
        output_text.insert(tk.END, result + '\n')
    
    output_text.insert(tk.END, "- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -\n", "bold")
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "Similarly to the basic Wolfrom part, pick a number for the pregearing and add it in the newly shown input field. Finally, click on the 'Summary' button to get a summary of the configurations you have chosen.", "bold")
    show_configuration_number_pre()

def summary():
    values = [entry.get() for entry in entry_boxes]
    hole_decision = hole_var.get()
    if hole_decision =='Yes':
        wolfrom_configuration_chosen = sorted_possible_solutions_wolfrom[int(values[20])-1]
        pregearing_configuration_chosen = sorted_final_solutions_pregearing[int(values[21])-1]
    else:
        wolfrom_configuration_chosen = sorted_possible_solutions_wolfrom[int(values[17])-1]
        pregearing_configuration_chosen = sorted_final_solutions_pregearing[int(values[18])-1]

    output_text.insert(tk.END, "Properties chosen configuration: ", "bold")
    output_text.insert(tk.END, "\n", "bold")
    output_text.insert(tk.END, "\n", "bold")
    Z_R2 = wolfrom_configuration_chosen[0]
    Z_R1 = wolfrom_configuration_chosen[1]
    Z_P1 = wolfrom_configuration_chosen[2]
    Z_P2 = wolfrom_configuration_chosen[3]
    GR_wolf = wolfrom_configuration_chosen[4]
    eff_wolf = wolfrom_configuration_chosen[5]
    eff_wolf_with_epsilon_alpha_equal_to_1 = wolfrom_configuration_chosen[6]
    b1 = wolfrom_configuration_chosen[7]
    b2 = wolfrom_configuration_chosen[8]
    Helix_angle_input_wolf = wolfrom_configuration_chosen[9]    
    Helix_angle_output_wolf = wolfrom_configuration_chosen[10]
    Epsilon_alpha_input_wolfrom_before_profile_shift = wolfrom_configuration_chosen[11]
    Epsilon_alpha_output_wolfrom_before_profile_shift = wolfrom_configuration_chosen[12]
    a_1 = wolfrom_configuration_chosen[13]
    a_2 = wolfrom_configuration_chosen[14]
        # Possible_solutions_wolfrom[i].append(1-f_R1[i])
        # Possible_solutions_wolfrom[i].append(1-f_R2[i])
        # Possible_solutions_wolfrom[i].append(20)
        # Possible_solutions_wolfrom[i].append(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_1[i]*np.pi/180))*180/np.pi)
        # Possible_solutions_wolfrom[i].append(np.arctan(np.tan(pressure_angle*np.pi/180)/np.cos(HA_2[i]*np.pi/180))*180/np.pi)
    F_t_R2 = wolfrom_configuration_chosen[20]
    F_t_R1 = wolfrom_configuration_chosen[21]
    F_t_P1 = wolfrom_configuration_chosen[22]
    F_t_P2 = wolfrom_configuration_chosen[23]
    SF_R2 = round(wolfrom_configuration_chosen[36],2)
    SF_R1 = round(wolfrom_configuration_chosen[37],2)
    SF_P1 = round(wolfrom_configuration_chosen[38],2)
    SF_P2 = round(wolfrom_configuration_chosen[39],2)
    SH_R2 = round(wolfrom_configuration_chosen[44],2)
    SH_R1 = round(wolfrom_configuration_chosen[45],2)
    SH_P1 = round(wolfrom_configuration_chosen[46],2)
    SH_P2 = round(wolfrom_configuration_chosen[47],2)

    Z_R_1 = pregearing_configuration_chosen[0]
    Z_P_1 = pregearing_configuration_chosen[1]  
    Z_s = pregearing_configuration_chosen[2]
    Z_P0 = pregearing_configuration_chosen[3]
    GR_pre = pregearing_configuration_chosen[4] 
    eff_pre = pregearing_configuration_chosen[5]
    eff_pre_with_epsilon_alpha_equal_to_1 = pregearing_configuration_chosen[6]
    b0 = pregearing_configuration_chosen[7] 
    b_1 = pregearing_configuration_chosen[8]
    Helix_angle_input_pre = pregearing_configuration_chosen[9]  
    Helix_angle_output_pre = pregearing_configuration_chosen[10]
    Epsilon_alpha_input_pre = pregearing_configuration_chosen[11]
    Epsilon_alpha_output_pre = pregearing_configuration_chosen[12]
    a_0 = pregearing_configuration_chosen[13]
    a__1 = pregearing_configuration_chosen[14]

    F_t_R_1 = pregearing_configuration_chosen[15]
    F_t_P_1 = pregearing_configuration_chosen[16]
    F_t_s = pregearing_configuration_chosen[17]
    F_t_P0 = pregearing_configuration_chosen[18]

    SF_R_1 = round(pregearing_configuration_chosen[27],2)
    SF_P_1 = round(pregearing_configuration_chosen[28],2)
    SF_S = round(pregearing_configuration_chosen[29],2)
    SF_P0 = round(pregearing_configuration_chosen[30],2)

    SH_R_1 = round(pregearing_configuration_chosen[39],2)
    SH_P_1 = round(pregearing_configuration_chosen[40],2)
    SH_S = round(pregearing_configuration_chosen[41],2)
    SH_P0 = round(pregearing_configuration_chosen[42],2)

    #after profile shift:
    x_R1 = round(wolfrom_configuration_chosen[52],4)
    x_P1 = round(wolfrom_configuration_chosen[53],4)
    a_1_PS = round(wolfrom_configuration_chosen[54],4)
    x_R2 = round(wolfrom_configuration_chosen[63],4)
    x_P2 = round(wolfrom_configuration_chosen[64],4)
    a_2_PS = round(wolfrom_configuration_chosen[65],4)
    SF_R1_PS = round(wolfrom_configuration_chosen[48],2)
    SF_R2_PS = round(wolfrom_configuration_chosen[59],2)
    SF_P1_PS = round(wolfrom_configuration_chosen[49],2)
    SF_P2_PS = round(wolfrom_configuration_chosen[60],2)

    SH_R1_PS = round(wolfrom_configuration_chosen[57],2)
    SH_R2_PS = round(wolfrom_configuration_chosen[70],2)
    SH_P1_PS = round(wolfrom_configuration_chosen[58],2)
    SH_P2_PS = round(wolfrom_configuration_chosen[71],2)

    x_s = round(pregearing_configuration_chosen[48],2)
    x_P0 = round(pregearing_configuration_chosen[47],2)
    x_P_1 = round(pregearing_configuration_chosen[58],2)
    x_R_1 = round(pregearing_configuration_chosen[59],2)
    a_0_PS = round(pregearing_configuration_chosen[49]*10**3,4)
    a__1_PS = round(pregearing_configuration_chosen[60]*10**3,4)
    SF_s_PS = round(pregearing_configuration_chosen[43],2)
    SF_P0_PS = round(pregearing_configuration_chosen[44],2)
    SF_P_1_PS = round(pregearing_configuration_chosen[55],2)
    SF_R_1_PS = round(pregearing_configuration_chosen[54],2)
    SH_s_PS = round(pregearing_configuration_chosen[52],2)
    SH_P0_PS = round(pregearing_configuration_chosen[53],2)
    SH_P_1_PS = round(pregearing_configuration_chosen[66],2)
    SH_R_1_PS = round(pregearing_configuration_chosen[65],2)

  
    output_text.insert(tk.END, "Before profile shift: ", "bold")
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "Results per gear: ", "bold")
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format(" ", "R2", "P2", "R1", "P1", "s", "P0", "R1'", "P1'")) 
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("NOT", Z_R2, Z_P2, Z_R1, Z_P1, Z_s, Z_P0, Z_R_1, Z_P_1))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("S_F", SF_R2, SF_P2, SF_R1, SF_P1, SF_S, SF_P0, SF_R_1, SF_P_1))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("S_H", SH_R2, SH_P2, SH_R1, SH_P1, SH_S, SH_P0, SH_R_1, SH_P_1))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_t [N]", round(F_t_R2, 2), round(F_t_P2, 2), round(F_t_R1, 2), round(F_t_P1, 2), round(F_t_s, 2), round(F_t_P0, 2), round(F_t_R_1, 2), round(F_t_P_1, 2)))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_a [N]", round(calculate_F_a(F_t_R2, Helix_angle_output_wolf), 2), round(calculate_F_a(F_t_P2, Helix_angle_output_wolf), 2), round(calculate_F_a(F_t_R1, Helix_angle_input_wolf), 2), round(calculate_F_a(F_t_P1, Helix_angle_input_wolf), 2), round(calculate_F_a(F_t_s, Helix_angle_input_pre), 2), round(calculate_F_a(F_t_P0, Helix_angle_input_pre), 2), round(calculate_F_a(F_t_R_1, Helix_angle_output_pre), 2), round(calculate_F_a(F_t_P_1, Helix_angle_output_pre), 2)))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_r [N]", round(calculate_F_r(F_t_R2, pressure_angle, Helix_angle_output_wolf), 2), round(calculate_F_r(F_t_P2, pressure_angle, Helix_angle_output_wolf), 2), round(calculate_F_r(F_t_R1, pressure_angle, Helix_angle_input_wolf), 2), round(calculate_F_r(F_t_P1, pressure_angle, Helix_angle_input_wolf), 2), round(calculate_F_r(F_t_s, pressure_angle, Helix_angle_input_pre), 2), round(calculate_F_r(F_t_P0, pressure_angle, Helix_angle_input_pre), 2), round(calculate_F_r(F_t_R_1, pressure_angle, Helix_angle_output_pre), 2), round(calculate_F_r(F_t_P_1, pressure_angle, Helix_angle_output_pre), 2)))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("F_n [N]", round(calculate_F_n(F_t_R2, pressure_angle), 2), round(calculate_F_n(F_t_P2, pressure_angle), 2), round(calculate_F_n(F_t_R1, pressure_angle), 2), round(calculate_F_n(F_t_P1, pressure_angle), 2), round(calculate_F_n(F_t_s, pressure_angle), 2), round(calculate_F_n(F_t_P0, pressure_angle), 2), round(calculate_F_n(F_t_R_1, pressure_angle), 2), round(calculate_F_n(F_t_P_1, pressure_angle), 2)))
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "Results per mesh: ", "bold")
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}\n".format(" ", "2", "1", "1'", "0")) 
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}\n".format("b [mm]", round(b2,2), round(b1,2), round(b_1,2), round(b0,2)))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}\n".format("a [mm]", round(a_2*10**3,2), round(a_1*10**3,2), round(a__1*10**3,2), round(a_0*10**3,2)))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}\n".format("Helix angle [°]", round(Helix_angle_output_wolf,2), round(Helix_angle_input_wolf,2), round(Helix_angle_output_pre,2), round(Helix_angle_input_pre,2)))
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "Results per part: ", "bold")
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "{:<30}{:<20}{:<20}{:<20}\n".format(" ", "Basic Wolfrom", "Pregearing", "Total'")) 
    output_text.insert(tk.END, "{:<30}{:<20}{:<20}{:<20}\n".format("GR", round(GR_wolf,2), round(GR_pre,2), round(GR_wolf*GR_pre,2)))
    output_text.insert(tk.END, "{:<30}{:<20}{:<20}{:<20}\n".format("Eff", round(eff_wolf,3), round(eff_pre,3), round(eff_wolf*eff_pre,3)))
    output_text.insert(tk.END, "{:<30}{:<20}{:<20}{:<20}\n".format("Eff if contact ratio = 1", round(eff_wolf_with_epsilon_alpha_equal_to_1,3), round(eff_pre_with_epsilon_alpha_equal_to_1,3), round(eff_pre_with_epsilon_alpha_equal_to_1*eff_wolf_with_epsilon_alpha_equal_to_1,3)))
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "After profile shift: ", "bold")
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format(" ", "R2", "P2", "R1", "P1", "s", "P0", "R1'", "P1'")) 
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("PS coeff", x_R2, x_P2, x_R1, x_P1, x_s, x_P0, x_R_1, x_P_1))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("S_F", SF_R2_PS, SF_P2_PS, SF_R1_PS, SF_P1_PS, SF_s_PS, SF_P0_PS, SF_R_1_PS, SF_P_1_PS))
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}{:<10}\n".format("S_H", SH_R2_PS, SH_P2_PS, SH_R1_PS, SH_P1_PS, SH_s_PS, SH_P0_PS, SH_R_1_PS, SH_P_1_PS))
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "Results per mesh: ", "bold")
    output_text.insert(tk.END, '\n')
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}\n".format(" ", "2", "1", "1'", "0")) 
    output_text.insert(tk.END, "{:<15}{:<10}{:<10}{:<10}{:<10}\n".format("a [mm]", a_2_PS*10**3, a_1_PS*10**3, a__1_PS, a_0_PS))

if __name__ == '__main__':
    create_gui()  # Run the GUI