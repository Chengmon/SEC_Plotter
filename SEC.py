import tkinter as tk
from tkinter import ttk, filedialog, messagebox
import matplotlib.pyplot as plt
import numpy as np
from pycorn import pc_uni6
import xml.etree.ElementTree as ET
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from ttkthemes import ThemedTk
from datetime import datetime

# Global variables
current_fig = None
coord_label = None
crosshair_vline = None
data_x1 = None
data_y1 = None
data_x2 = None
data_y2 = None
uv1_name = "UV1"
uv2_name = "UV2"
frac_x = None
frac_labels = None
creation_date = "N/A"  

# Function to update coordinates, crosshair, and fraction on hover
def update_coords(event):
    global coord_label, crosshair_vline, data_x1, data_y1, data_x2, data_y2, uv1_name, uv2_name, frac_x, frac_labels
    if event.inaxes and current_fig and data_x1 is not None and data_x2 is not None and frac_x is not None:
        x = event.xdata
        # Find the closest x-value in data_x1 and data_x2 for UV values
        idx1 = np.argmin(np.abs(data_x1 - x))
        idx2 = np.argmin(np.abs(data_x2 - x))
        nearest_x1 = data_x1[idx1]
        nearest_y1 = data_y1[idx1]
        nearest_x2 = data_x2[idx2]
        nearest_y2 = data_y2[idx2]
        
        # Use the x-value closest to the cursor for the crosshair
        nearest_x = nearest_x1 if abs(nearest_x1 - x) < abs(nearest_x2 - x) else nearest_x2
        
        # Find the fraction range that x falls into
        nearest_frac = "N/A"
        for i in range(len(frac_x)):
            start = frac_x[i]
            end = frac_x[i + 1] if i + 1 < len(frac_x) else float('inf')
            if start <= x < end:
                nearest_frac = frac_labels[i]
                break
        
        # Update coordinates label with UV names and fraction
        coord_label.config(text=f"Elution: {nearest_x:.2f} mL, UV{uv1_name}: {nearest_y1:.2f} mAU, UV{uv2_name}: {nearest_y2:.2f} mAU, Fraction: {nearest_frac}")
        
        # Update vertical crosshair position
        if crosshair_vline is None:
            crosshair_vline = event.inaxes.axvline(nearest_x, color='gray', linestyle='--', linewidth=0.5)
        else:
            crosshair_vline.set_xdata([nearest_x, nearest_x])
        canvas.draw()
    else:
        coord_label.config(text=f"Elution: N/A, UV{uv1_name}: N/A, UV{uv2_name}: N/A, Fraction: N/A")
        if crosshair_vline is not None:
            crosshair_vline.set_xdata([np.nan, np.nan])
            canvas.draw()

# Function to process and plot the data
def process_and_plot(file_path, xlim_left, xlim_right, ylim_bottom, ylim_top, canvas):
    global current_fig, data_x1, data_y1, data_x2, data_y2, crosshair_vline, uv1_name, uv2_name, frac_x, frac_labels, creation_date, date_label
    try:
        my_res_file = pc_uni6(file_path)
        my_res_file.load()
        my_res_file.xml_parse()

        x1 = np.array(my_res_file['Chrom.1_1_True']['CoordinateData.Volumes'])
        y1 = np.array(my_res_file['Chrom.1_1_True']['CoordinateData.Amplitudes'])
        x2 = np.array(my_res_file['Chrom.1_2_True']['CoordinateData.Volumes'])
        y2 = np.array(my_res_file['Chrom.1_2_True']['CoordinateData.Amplitudes'])
        inject1 = my_res_file['Injection']['data'][0][0]

        # Store data globally for coordinate lookup
        data_x1 = x1 - inject1
        data_y1 = y1
        data_x2 = x2 - inject1
        data_y2 = y2

        x_data = np.concatenate([x1, x2])
        y_data = np.concatenate([y1, y2])
        xlim_left_default = float(np.min(x_data))
        xlim_right_default = float(np.max(x_data))
        ylim_bottom_default = float(np.min(y_data))
        ylim_top_default = float(np.max(y_data)*1.1)

        if not entry_xlim_left.get():
            entry_xlim_left.delete(0, tk.END)
            entry_xlim_left.insert(0, f"{xlim_left_default:.2f}")
        if not entry_xlim_right.get():
            entry_xlim_right.delete(0, tk.END)
            entry_xlim_right.insert(0, f"{xlim_right_default:.2f}")
        if not entry_ylim_bottom.get():
            entry_ylim_bottom.delete(0, tk.END)
            entry_ylim_bottom.insert(0, f"{ylim_bottom_default:.2f}")
        if not entry_ylim_top.get():
            entry_ylim_top.delete(0, tk.END)
            entry_ylim_top.insert(0, f"{ylim_top_default:.2f}")

        # Parse UV names and creation date from XML
        chrom = ET.fromstring(my_res_file['Result.xml'])
        uv1_name = chrom.find('.//*[Keyword1="UV1"]/Keyword2').text
        uv2_name = chrom.find('.//*[Keyword1="UV2"]/Keyword2').text
        column = chrom.find('.//*[Keyword1="Column type"]/Keyword2').text
        creation_raw = chrom.find('./Created').text  # e.g., "2024-01-31T13:21:50.208"
        creation_date = datetime.strptime(creation_raw.split('.')[0], "%Y-%m-%dT%H:%M:%S").strftime("%d-%b-%Y %H:%M")


        # Update the date label
        date_label.config(text=f"Purification date: {creation_date}")

        # Extract fraction data from my_res_file
        frac_data = my_res_file['Fractions']['data']
        frac_x = np.array([x for x, _ in frac_data])
        frac_labels = [label for _, label in frac_data]

        fig = plt.Figure(dpi=110)
        ax = fig.add_subplot(111)
        ax.plot(data_x1, y1, 'b-', linewidth=0.8, label=f'UV {uv1_name}')
        ax.plot(data_x2, y2, 'r-', linewidth=0.8, label=f'UV {uv2_name}')
        ax.set_xlabel('Elution volume [mL]')
        ax.set_ylabel('Absorbance [mAu]')


        try:
            ax.set_xlim(left=float(xlim_left) if xlim_left else xlim_left_default, 
                       right=float(xlim_right) if xlim_right else xlim_right_default)
            ax.set_ylim(bottom=float(ylim_bottom) if ylim_bottom else ylim_bottom_default, 
                       top=float(ylim_top) if ylim_top else ylim_top_default)
        except ValueError:
            messagebox.showwarning("Warning", "Invalid limit values; using data range.")
            ax.set_xlim(left=xlim_left_default, right=xlim_right_default)
            ax.set_ylim(bottom=ylim_bottom_default, top=ylim_top_default)

        ax.set_title(column, loc='left', fontsize=7)
        fig.suptitle('Column', horizontalalignment = 'left', x=fig.subplotpars.left, y=0.96, fontweight=1000)
        ax.legend(loc='upper right', bbox_to_anchor=(1, 1.11), 
                  fontsize='x-small', frameon=False)
        ax.minorticks_on()

        # Reset crosshair
        crosshair_vline = None

        current_fig = fig
        canvas.figure = fig
        canvas.draw()

        # Connect the motion event to update coordinates and crosshair
        canvas.mpl_connect('motion_notify_event', update_coords)

    except Exception as e:
        messagebox.showerror("Error", f"Failed to process file: {str(e)}")

# Function to save the plot
def save_plot():
    global current_fig
    if current_fig is None:
        messagebox.showwarning("Warning", "No plot to save. Generate a plot first!")
        return
    save_path = filedialog.asksaveasfilename(
        title="Save Plot As",
        defaultextension=".png",
        filetypes=(("PNG files", "*.png"), ("JPEG files", "*.jpg"), ("All files", "*.*"))
    )
    if save_path:
        try:
            current_fig.savefig(save_path, dpi=300, bbox_inches='tight')
            messagebox.showinfo("Success", f"Plot saved to {save_path}")
        except Exception as e:
            messagebox.showerror("Error", f"Failed to save plot: {str(e)}")

# Function to browse and select a file
def browse_file():
    file_path = filedialog.askopenfilename(
        title="Select Chromatography ZIP File",
        filetypes=(("ZIP files", "*.zip"), ("All files", "*.*"))
    )
    if file_path:
        selected_file.set(file_path)
        file_label.config(text=f"Selected: {file_path}")
        entry_xlim_left.delete(0, tk.END)
        entry_xlim_right.delete(0, tk.END)
        entry_ylim_bottom.delete(0, tk.END)
        entry_ylim_top.delete(0, tk.END)
    else:
        file_label.config(text="No file selected")

# Function to plot with user-defined limits
def plot_file():
    if selected_file.get():
        xlim_left = entry_xlim_left.get()
        xlim_right = entry_xlim_right.get()
        ylim_bottom = entry_ylim_bottom.get()
        ylim_top = entry_ylim_top.get()
        process_and_plot(selected_file.get(), xlim_left, xlim_right, 
                         ylim_bottom, ylim_top, canvas)
    else:
        messagebox.showwarning("Warning", "Please select a file first!")

# Create the GUI with ttkthemes
window = ThemedTk(theme="arc")
window.title("AKTA Pure SEC Plotter")
window.geometry("1000x600")

# Variable to store the selected file path
selected_file = tk.StringVar()

# Left frame for controls
left_frame = ttk.Frame(window, padding=10, width=300)
left_frame.pack(side=tk.LEFT, fill=tk.Y)

# File selection section
file_label = ttk.Label(left_frame, text="No file selected", wraplength=250, font=("Segoe UI", 10))
file_label.pack(pady=5)

browse_button = ttk.Button(left_frame, text="Select ZIP File", command=browse_file, style="Custom.TButton")
browse_button.pack(pady=5)

# Plot limits section
limits_frame = ttk.LabelFrame(left_frame, text="Plot Limits", padding=10)
limits_frame.pack(fill=tk.X, pady=10)

ttk.Label(limits_frame, text="X-Limit Left:", font=("Segoe UI", 9)).grid(row=0, column=0, padx=5, pady=2)
entry_xlim_left = ttk.Entry(limits_frame, width=10, font=("Segoe UI", 9))
entry_xlim_left.grid(row=0, column=1, padx=5, pady=2)

ttk.Label(limits_frame, text="X-Limit Right:", font=("Segoe UI", 9)).grid(row=1, column=0, padx=5, pady=2)
entry_xlim_right = ttk.Entry(limits_frame, width=10, font=("Segoe UI", 9))
entry_xlim_right.grid(row=1, column=1, padx=5, pady=2)

ttk.Label(limits_frame, text="Y-Limit Bottom:", font=("Segoe UI", 9)).grid(row=2, column=0, padx=5, pady=2)
entry_ylim_bottom = ttk.Entry(limits_frame, width=10, font=("Segoe UI", 9))
entry_ylim_bottom.grid(row=2, column=1, padx=5, pady=2)

ttk.Label(limits_frame, text="Y-Limit Top:", font=("Segoe UI", 9)).grid(row=3, column=0, padx=5, pady=2)
entry_ylim_top = ttk.Entry(limits_frame, width=10, font=("Segoe UI", 9))
entry_ylim_top.grid(row=3, column=1, padx=5, pady=2)

# Button section
button_frame = ttk.Frame(left_frame)
button_frame.pack(pady=10)

plot_button = ttk.Button(button_frame, text="Generate Plot", command=plot_file, style="Custom.TButton")
plot_button.pack(pady=5)

save_button = ttk.Button(button_frame, text="Save Plot", command=save_plot, style="Custom.TButton")
save_button.pack(pady=5)

# Date label at the bottom of left frame
date_label = ttk.Label(left_frame, text=f"Purification date: {creation_date}", font=("Segoe UI", 9))
date_label.pack(side=tk.BOTTOM, pady=10, anchor="sw")

# Right frame for plot and coordinates
right_frame = ttk.Frame(window, relief="sunken", borderwidth=1)
right_frame.pack(side=tk.RIGHT, fill=tk.BOTH, expand=True, padx=10, pady=10)

canvas = FigureCanvasTkAgg(plt.Figure(figsize=(6, 4), dpi=100), master=right_frame)
canvas.get_tk_widget().pack(fill=tk.BOTH, expand=True)

# Coordinates label below the plot
coord_label = ttk.Label(right_frame, text=f"Elution: N/A, UV{uv1_name}: N/A, UV{uv2_name}: N/A, Fraction: N/A", font=("Segoe UI", 9))
coord_label.pack(pady=5)

# Custom button style
style = ttk.Style()
style.configure("Custom.TButton", font=("Segoe UI", 10, "bold"), padding=6, background="#4CAF50", foreground="black")
style.map("Custom.TButton", background=[("active", "#45a049")]) 

# Start the GUI
window.mainloop()
