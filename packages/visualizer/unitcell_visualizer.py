""" visualize the unit cell as done in VESTA
""" 

import sys
import os

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import Dans_Diffraction as dif
import matplotlib.cm as cm
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt

# Import Dans_Diffraction internal functions for plotting
from Dans_Diffraction import functions_general as fg
from Dans_Diffraction import functions_plotting as fp
from Dans_Diffraction import functions_crystallography as fc

# Apply the colormap patch for Dans_Diffraction compatibility
_old_ensure_cmap = cm._ensure_cmap

def _ensure_cmap_patched(cmap):
    if isinstance(cmap, np.ndarray):
        return mcolors.ListedColormap(cmap)
    return _old_ensure_cmap(cmap)

cm._ensure_cmap = _ensure_cmap_patched


class UnitcellVisualizer(FigureCanvas):
    def __init__(self, width=4, height=4, dpi=100):
        self.fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = self.fig.add_subplot(111, projection='3d')
        self.cif_file_path = None
        self._is_initialized = False
        self.crystal = None
        super().__init__(self.fig)

        # Set background color to white
        self.fig.patch.set_facecolor("white")
        self.axes.set_facecolor("white")

        # Set initial view
        self.axes.view_init(elev=39, azim=75)

        # Clean axes like Dans_Diffraction
        self.axes.set_axis_off()

    @property
    def cif_file_path(self):
        return self._cif_file_path

    @cif_file_path.setter
    def cif_file_path(self, cif_file_path: str):
        # check if the file exists (allow None for initialization)
        if cif_file_path is not None and not os.path.exists(cif_file_path):
            raise FileNotFoundError(f"File {cif_file_path} does not exist.")
        self._cif_file_path = cif_file_path

    @property
    def is_initialized(self):
        return self._is_initialized

    @is_initialized.setter
    def is_initialized(self, is_initialized: bool):
        self._is_initialized = is_initialized

    
    def set_parameters(self, params: dict):
        cif_file_path = params["cif_file"]
        self.cif_file_path = cif_file_path
        try:
            self.crystal = dif.Crystal(self.cif_file_path)
            self.is_initialized = True
        except Exception as e:
            print(f"Error loading crystal structure from {cif_file_path}: {e}")
            self.is_initialized = False

    def visualize_unitcell(self, show_labels=False):
        """ read the atom position from the cif file and visualize the unit cell using integrated Dans_Diffraction logic
        """ 
        if not self.is_initialized or self.crystal is None:
            print("Unit cell visualizer is not initialized. Please set parameters first.")
            return
        
        # Clear the current plot
        self.axes.clear()
        
        # Set background color to white and clean axes like Dans_Diffraction
        self.axes.set_facecolor("white")
        self.axes.set_axis_off()  # Ensure clean look after clearing
        
        try:
            # Direct implementation of Dans_Diffraction plot_crystal logic
            self._plot_crystal_direct(show_labels=show_labels)
            
            # Set title with crystal info
            #title = f"Unit Cell: {os.path.basename(self.cif_file_path)}"
            #if hasattr(self.crystal, 'name') and self.crystal.name:
            #    title = f"Unit Cell: {self.crystal.name}"
            #title = f"Unit Cell"
            #self.axes.set_title(title, fontsize=12, pad=20)
            
        except Exception as e:
            print(f"Error visualizing unit cell: {e}")
            # Fallback: show error message on plot
            self.axes.text(0.5, 0.5, 0.5, f"Error: {str(e)}", 
                          transform=self.axes.transAxes, 
                          ha='center', va='center')
        
        # Set initial view
        self.axes.view_init(elev=30, azim=55)
        
        # Refresh the canvas
        self.draw()
    
    def _plot_crystal_direct(self, show_labels=False):
        """
        Direct implementation of Dans_Diffraction plot_crystal method for our Qt axes
        Based on Dans_Diffraction.classes_plotting.Plotting.plot_crystal
        """
        # Generate lattice
        tol = 0.05
        uvw, element, label, occ, uiso, mxmymz = self.crystal.Structure.generate_lattice(1, 1, 1)
        
        # Split atom types, color & radii
        labels, idx, invidx = np.unique(label, return_index=True, return_inverse=True)
        types = element[idx]
        colors = plt.cm.rainbow(np.linspace(0, 1, len(types)))
        sizes = fc.atom_properties(types, 'Radii')
        
        # Get atomic positions
        R = self.crystal.Cell.calculateR(uvw)
        cen = np.mean(R, axis=0)
        R = R - cen
        I = np.all(np.hstack([uvw<(1+tol), uvw>(0-tol), occ.reshape([-1,1])>0.2]), 1)
        
        # Magnetic vectors
        V = self.crystal.Cell.calculateR(mxmymz/np.asarray(self.crystal.Cell.lp()[:3]))
        # Get axis limits based on lattice (needed for coordinate arrows positioning)
        lim = np.max(self.crystal.Cell.lp()[:3])
        # Loop over each atom type and plot them
        for n in range(len(types)):
            # don't plot unoccupied positions
            tot_occ = np.array([occ[m] for m in range(len(R)) if invidx[m] == n])
            if sum(tot_occ) == 0: 
                continue

            xyz = np.array([R[m, :] for m in range(len(R)) if invidx[m] == n])
            iii = np.array([I[m] for m in range(len(R)) if invidx[m] == n])
            col = np.tile(colors[n], (len(xyz[iii, :]), 1))
            
            # Plot atoms as scatter points
            self.axes.scatter(xyz[iii, 0], xyz[iii, 1], xyz[iii, 2], 
                            s=1.5 * sizes[n], c=col, label=labels[n], alpha=0.85,
                            edgecolors='white', linewidth=0.1)
            
            # Plot magnetic vectors (arrows) if they exist
            for m in range(len(R)): 
                if invidx[m] == n and I[m]:
                    xyz_pos = R[m, :]
                    vec = V[m, :]
                    if fg.mag(vec) < 0.1: 
                        continue
                    # Simple arrow implementation (could be enhanced)
                    start = xyz_pos - vec / 2
                    end = xyz_pos + vec / 2
                    self.axes.plot([start[0], end[0]], [start[1], end[1]], [start[2], end[2]], 
                                 'r-', linewidth=3, alpha=0.8)
        
        # Labels
        if show_labels:
            uvw_st, type_st, label_st, occ_st, uiso_st, mxmymz_st = self.crystal.Structure.get()
            R_st = self.crystal.Cell.calculateR(uvw_st) - cen
            for n in range(len(R_st)):
                self.axes.text(R_st[n, 0], R_st[n, 1], R_st[n, 2], 
                             '%2d: %s' % (n, label_st[n]), fontsize=8,
                             bbox=dict(boxstyle="round,pad=0.2", facecolor='white', alpha=0.8))
        
        # Create cell box
        uvw_box = np.array([[0., 0, 0], [1, 0, 0], [1, 0, 1], [1, 1, 1], [1, 1, 0], [0, 1, 0], [0, 1, 1],
                           [0, 0, 1], [1, 0, 1], [1, 0, 0], [1, 1, 0], [1, 1, 1], [0, 1, 1], [0, 1, 0], [0, 0, 0],
                           [0, 0, 1]])
        bpos = self.crystal.Cell.calculateR(uvw_box) - cen
        self.axes.plot(bpos[:, 0], bpos[:, 1], bpos[:, 2], c='gray', linewidth=1, alpha=0.85)  # cell box
        
        # Plot lattice vector arrows (a=red, b=green, c=blue)
        # a vector (red)
        self.axes.plot([bpos[0, 0], bpos[1, 0]], [bpos[0, 1], bpos[1, 1]], [bpos[0, 2], bpos[1, 2]], 
                      'r-', linewidth=4, alpha=0.0)
        # b vector (green)  
        self.axes.plot([bpos[0, 0], bpos[5, 0]], [bpos[0, 1], bpos[5, 1]], [bpos[0, 2], bpos[5, 2]], 
                      'g-', linewidth=4, alpha=0.0)
        # c vector (blue)
        self.axes.plot([bpos[0, 0], bpos[7, 0]], [bpos[0, 1], bpos[7, 1]], [bpos[0, 2], bpos[7, 2]], 
                      'b-', linewidth=4, alpha=0.0)
        
        
        
        # Create coordinate basis arrows outside the cell box
        # Position them in a corner away from the main structure
        arrow_origin = np.array([-lim*0.4, -lim*0.4, -lim*0.4])  # Bottom-left-front corner
        arrow_scale = lim * 0.15  # Scale factor for arrow length
        
        # Get normalized lattice vectors from the crystal
        a_vec = np.array([1, 0, 0])
        b_vec = np.array([0, 1, 0]) 
        c_vec = np.array([0, 0, 1])
        
        # Transform to real space coordinates
        a_real = self.crystal.Cell.calculateR(a_vec.reshape(1, -1))[0]
        b_real = self.crystal.Cell.calculateR(b_vec.reshape(1, -1))[0]
        c_real = self.crystal.Cell.calculateR(c_vec.reshape(1, -1))[0]
        
        # Normalize and scale the vectors
        a_norm = a_real / np.linalg.norm(a_real) * arrow_scale
        b_norm = b_real / np.linalg.norm(b_real) * arrow_scale
        c_norm = c_real / np.linalg.norm(c_real) * arrow_scale
        
        # Draw coordinate arrows with labels
        # a vector (red)
        self.axes.quiver(arrow_origin[0], arrow_origin[1], arrow_origin[2]+ 2*bpos[7, 2],
                        a_norm[0], a_norm[1], a_norm[2],
                        color='dodgerblue', arrow_length_ratio=0.15, linewidth=2, alpha=0.9)
        self.axes.text(arrow_origin[0] + a_norm[0]*1.2, arrow_origin[1] + a_norm[1]*1, 
                      arrow_origin[2] + a_norm[2]*1.2 + 2*bpos[7, 2], 'a', color='dodgerblue', fontsize=14, fontweight='bold')
        
        # b vector (green)
        self.axes.quiver(arrow_origin[0], arrow_origin[1], arrow_origin[2]+ 2*bpos[7, 2],
                        b_norm[0], b_norm[1], b_norm[2],
                        color='dodgerblue', arrow_length_ratio=0.15, linewidth=2, alpha=0.9)
        self.axes.text(arrow_origin[0] + b_norm[0]*1.2, arrow_origin[1] + b_norm[1]*1.2, 
                      arrow_origin[2] + b_norm[2]*1.2 + 2*bpos[7, 2], 'b', color='dodgerblue', fontsize=14, fontweight='bold')
        
        # c vector (blue)
        self.axes.quiver(arrow_origin[0], arrow_origin[1], arrow_origin[2]+ 2*bpos[7, 2],
                        c_norm[0], c_norm[1], -c_norm[2],
                        color='dodgerblue', arrow_length_ratio=0.15, linewidth=2, alpha=0.9)
        self.axes.text(arrow_origin[0] + c_norm[0]*1.2, arrow_origin[1] + c_norm[1]*1.2, 
                      arrow_origin[2] - c_norm[2]*1.2 + 2*bpos[7, 2], 'c', color='dodgerblue', fontsize=14, fontweight='bold')
        
        # Set axis limits based on lattice
        self.axes.set_xlim(-lim/2, lim/2)
        self.axes.set_ylim(-lim/2, lim/2)
        self.axes.set_zlim(-lim/2, lim/2)
        
        # Turn off axis (removes background, grids, ticks, labels) - like Dans_Diffraction
        self.axes.set_axis_off()
        # change view angles
        # Add legend positioned closer to the unit cell
        self.axes.legend(fontsize=13, frameon=False, loc='upper right', 
                        bbox_to_anchor=(0.98, 0.98), handletextpad=0.3, 
                        handlelength=1.5, columnspacing=0.5, 
                        fancybox=False, shadow=False)
        self.fig.tight_layout()

if __name__ == "__main__":
    import matplotlib
    import matplotlib.pyplot as plt
    
    # Set matplotlib to use a non-interactive backend for testing
    matplotlib.use('Agg')  # Use Agg backend for saving files
    
    print(f"Using matplotlib backend: {matplotlib.get_backend()}")
    print("Testing UnitcellVisualizer...")
    
    # Create visualizer
    viz = UnitcellVisualizer(width=10, height=8)
    
    # Test with both available CIF files  
    cif_files = ['data/nacl.cif', 'data/LSCO.cif']
    
    for cif_file in cif_files:
        if os.path.exists(cif_file):
            print(f"Processing {cif_file}...")
            
            # Set CIF file path
            viz.set_parameters({'cif_file': cif_file})
            
            # Create visualization
            viz.visualize_unitcell()
            
            # Create output filename
            base_name = os.path.splitext(os.path.basename(cif_file))[0]
            output_file = f'figures/unitcell_{base_name}.png'
            
            # Ensure figures directory exists
            os.makedirs('figures', exist_ok=True)
            
            # Save the plot to file
            viz.fig.savefig(output_file, dpi=150, bbox_inches='tight', 
                           facecolor='white', edgecolor='none')
            print(f"Saved unit cell visualization to {output_file}")
            
        else:
            print(f"CIF file not found: {cif_file}")
    
    print("Test completed successfully!")
