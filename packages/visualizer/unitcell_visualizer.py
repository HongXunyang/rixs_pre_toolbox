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

from ..utils import angle_to_matrix, get_rotation

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

    def visualize_unitcell(self, show_labels=False, is_clear=True):
        """ read the atom position from the cif file and visualize the unit cell using integrated Dans_Diffraction logic
        """ 
        if not self.is_initialized or self.crystal is None:
            print("Unitcell visualizer not initialized.")
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
                            s=1 * sizes[n], c=col, label=labels[n], alpha=1,
                            edgecolors='white', linewidth=0.1)
           
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
        self.axes.quiver(arrow_origin[0], arrow_origin[1], arrow_origin[2]+ 1.6*bpos[7, 2],
                        a_norm[0], a_norm[1], a_norm[2],
                        color='dodgerblue', arrow_length_ratio=0.15, linewidth=2, alpha=0.9)
        self.axes.text(arrow_origin[0] + a_norm[0]*1.2, arrow_origin[1] + a_norm[1]*1, 
                      arrow_origin[2] + a_norm[2]*1.2 + 1.6*bpos[7, 2], 'a', color='dodgerblue', fontsize=14, fontweight='bold')
        
        self.axes.quiver(arrow_origin[0], arrow_origin[1], arrow_origin[2]+ 1.6*bpos[7, 2],
                        b_norm[0], b_norm[1], b_norm[2],
                        color='dodgerblue', arrow_length_ratio=0.15, linewidth=2, alpha=0.9)
        self.axes.text(arrow_origin[0] + b_norm[0]*1.2, arrow_origin[1] + b_norm[1]*1.2, 
                      arrow_origin[2] + b_norm[2]*1.2 + 1.6*bpos[7, 2], 'b', color='dodgerblue', fontsize=14, fontweight='bold')
        
        self.axes.quiver(arrow_origin[0], arrow_origin[1], arrow_origin[2]+ 1.6*bpos[7, 2],
                        c_norm[0], c_norm[1], c_norm[2],
                        color='dodgerblue', arrow_length_ratio=0.15, linewidth=2, alpha=0.9)
        self.axes.text(arrow_origin[0] + c_norm[0]*1.2, arrow_origin[1] + c_norm[1]*1.2, 
                      arrow_origin[2] + c_norm[2]*1.2 + 1.6*bpos[7, 2], 'c', color='dodgerblue', fontsize=14, fontweight='bold')
        
        # Set axis limits based on lattice
        self.axes.set_xlim(-lim/2, lim/2)
        self.axes.set_ylim(-lim/2, lim/2)
        self.axes.set_zlim(-lim/2, lim/2)
        
        # Turn off axis (removes background, grids, ticks, labels) - like Dans_Diffraction
        #self.axes.set_axis_off()
        # change view angles
        # Add legend positioned closer to the unit cell
        self.axes.legend(fontsize=10, frameon=False, loc='upper left', 
                        bbox_to_anchor=(0.98, 0.98), handletextpad=0.3, 
                        handlelength=1.5, columnspacing=0.1, 
                        fancybox=False, shadow=False)
        self.fig.tight_layout()




    def visualize_scattering_geometry(self, scattering_angles=None, is_clear=False):
        """ plotting the scattering plane and the beam as done in the ScatteringVisualizer
        """ 
        if not self.is_initialized or self.crystal is None:
            print("Unitcell visualizer not initialized.")
            return 
        if is_clear:
            # Clear previous plot
            self.axes.clear()
        # Plot the x-ray beam
        if scattering_angles is None:
            scattering_angles = {
                "theta": 50,
                "tth": 150,
                "phi": 0,
                "chi": 0,
            }
        tth, theta, phi, chi = scattering_angles.get("tth", 150), scattering_angles.get("theta", 50), scattering_angles.get("phi", 0), scattering_angles.get("chi", 0)
        a_vec = np.array([1, 0, 0])
        b_vec = np.array([0, 1, 0]) 
        c_vec = np.array([0, 0, 1])
        
        # Transform to real space coordinates
        a_real = self.crystal.Cell.calculateR(a_vec.reshape(1, -1))[0]
        b_real = self.crystal.Cell.calculateR(b_vec.reshape(1, -1))[0]
        c_real = self.crystal.Cell.calculateR(c_vec.reshape(1, -1))[0]



        # Get scale factor from unit cell if available
        if self.is_initialized and self.crystal is not None:
            lim = np.max(self.crystal.Cell.lp()[:3])
            scale_factor = lim / 2 * 1.5  # Match the unit cell scale
        else:
            scale_factor = 1.5  # Default scale if no crystal loaded
        # Plot the scattering plane - scale it to match unit cell, adjusted for beam from -y
        plane_width = 0.25 * scale_factor  # narrower in x (rotated 90°)
        plane_height_bottom = -0.75 * scale_factor  # extends more in y direction
        plane_height_top = 1.25 * scale_factor
        
        x_basis = a_real / np.linalg.norm(a_real)
        # y_basis = cross(a,b)
        z_basis = np.cross(a_real, b_real) / np.linalg.norm(np.cross(a_real, b_real))
        y_basis = np.cross(z_basis, x_basis) / np.linalg.norm(np.cross(z_basis, x_basis))

        # Define vertices in local plane coordinates (x-y plane, z=0)
        local_vertices = np.array([
            [plane_width*4, plane_height_bottom, 0],   # bottom right in local coords
            [-plane_width, plane_height_bottom, 0],  # bottom left in local coords
            [plane_width*4, plane_height_top, 0],      # top right in local coords
            [-plane_width, plane_height_top, 0],     # top left in local coords
        ])
        
        # Transform vertices to real space using crystal basis vectors
        scatter_plane_vertices = np.zeros_like(local_vertices)
        for i, vertex in enumerate(local_vertices):
            # Transform each vertex using the basis vectors
            # x component uses x_basis, y component uses y_basis, z component uses z_basis (normal to plane)
            scatter_plane_vertices[i] = (vertex[0] * x_basis + 
                                       vertex[1] * y_basis + 
                                       vertex[2] * z_basis)
        # coordinate change matrix: change to scattering plane system
        ccm = np.array([x_basis, y_basis, z_basis]).T
        scatter_plane_vertices = _rotate_vertices_wrt_plane(scatter_plane_vertices, ccm, phi, chi)
        scatter_plane_faces = np.array([[0, 1, 3, 2]])  # single face
        self.axes.add_collection3d(
            Poly3DCollection(
                scatter_plane_vertices[scatter_plane_faces],
                facecolors=[0.3510, 0.7850, 0.9330],  # light blue
                edgecolors=[0.7, 0.7, 0.7],
                alpha=0.3,
            )
        )



        # Plot incident beam (k_in) - coming from -y direction in local coords
        offset = 0
        k_in_length = 1.3 * scale_factor
        # In local scattering plane coords: -y_basis direction (rotated 90° from -x)
        k_in_vec = - k_in_length * y_basis
        k_in_vec = rotate_vector(k_in_vec, ccm, theta, phi, chi)
        # Draw colored arrow on top
        self.axes.quiver(
            -k_in_vec[0],
            -k_in_vec[1] - offset,
            -k_in_vec[2],
            k_in_vec[0],
            k_in_vec[1],
            k_in_vec[2],
            color=(191 / 255, 44 / 255, 0),
            alpha=1,
            linewidth=5,
            arrow_length_ratio=0.2,
            zorder=10,
        )

        # Plot scattered beam (k_out) - in x-y plane, rotated 90° to match -y incident beam
        k_out_length = 1.3 * scale_factor
        # Rotated 90° in local coords: incident from -y, scattered at angle tth
        k_out_vec = ccm @ np.array([np.sin(np.radians(tth)), -np.cos(np.radians(tth)), 0]) * k_out_length
        k_out_vec = rotate_vector(k_out_vec, ccm, theta, phi, chi)
        # Draw colored arrow on top
        self.axes.quiver(
            0,
            0 + offset,
            0,
            k_out_vec[0],
            k_out_vec[1],
            k_out_vec[2],
            color=(2 / 255, 78 / 255, 191 / 255),
            linewidth=5,
            arrow_length_ratio=0.2,
            zorder=10,
        )

        # Set axis limits to match unit cell scale
        if self.is_initialized and self.crystal is not None:
            self.axes.set_xlim(-lim/2, lim/2)
            self.axes.set_ylim(-lim/2, lim/2)
            self.axes.set_zlim(-lim/2, lim/2)
        else:
            # Fallback to default limits
            self.axes.set_xlim(-1, 1)
            self.axes.set_ylim(-1, 1)
            self.axes.set_zlim(-1, 1)

        #self.axes.set_axis_off()
        # Update the canvas
        self.draw()      
    
    def clear_plot(self):
        """ clear the plot """
        self.axes.clear()
        self.draw()
    

def rotate_vector(vec, ccm, theta, phi, chi):
    """Convert angles theta, phi, chi to rotation matrix. This rotate the beam or the scattering
    plane, not the sample. 
    Pay attention to the direction of the rotation.
    
    Updated for x-y scattering plane (z-axis is normal to scattering plane).

    Args:
        vec (np.ndarray): vector to rotate
        ccm (np.ndarray): coordinate change matrix
        theta (float): rotation about the z-axis in degrees, right-hand rule
        phi (float): rotation about the x-axis in degrees, right-hand rule
        chi (float): rotation about the y-axis in degrees, right-hand rule

    Returns:
        rotated_vec (np.ndarray): Rotated vector
    """

    matrix = angle_to_matrix(theta, phi, chi)
    matrix = ccm @ matrix.T @ ccm.T
    return matrix @ vec

def _rotate_vertices_wrt_plane(vertices, ccm, phi, chi):
    """Rotate the vertices of the sample with respect to the scattering plane
    
    Updated for x-y scattering plane (z-axis is normal to scattering plane).
    phi: rotation about y-axis
    chi: rotation about x-axis
    """

    rotation_matrix = get_rotation(phi, chi)
    rotation_matrix = ccm @ rotation_matrix.T @ ccm.T
    vertices = np.array(vertices)
    for i, vertex in enumerate(vertices):
        vertices[i] = rotation_matrix @ vertex
    return vertices






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


