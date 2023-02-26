#Required imports
import sys
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
import tkinter.font as tkfont

#Default arg values
arg_entries = False
arg_split = False
arg_textonly = False

#Read command line args
for arg in sys.argv[1:]:
    if not arg.startswith('-'):
        continue
    if 'e' in arg:
        arg_entries = True
    if 's' in arg:
        arg_split = True
    if 't' in arg:
        arg_textonly = True

#Try to import FigureCanvasTkAgg to composite pyplot into tk GUI window
if not arg_textonly:
    try:
        from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
    except ImportError as err:
        print('Error importing FigureCanvasTkAgg - Using simple GUI instead')
        arg_textonly = True

#Have to split windows and use entries if using text-only GUI
if arg_textonly:
    arg_split = True
    arg_entries = True

#If splitting windows, need to modify global setting in matplotlib
if arg_split:
    import matplotlib as mpl
    mpl.use('TkAgg')
    mpl.rcParams['toolbar'] = 'None'

#Try to import matplotlib's RangeSlider widget if needed
if not arg_entries:
    try:
        from matplotlib.widgets import RangeSlider
    except ImportError as err:
        print('Error importing RangeSlider - Using Tkinter entry widgets instead')
        arg_entries = True

#Aspect ratio 2/sqrt(3)
ratio = 2 / 3 ** 0.5

class Snowflake:
    #Default properties
    min_size = 7
    bounds = [None] * 4
    branch = False
    zooming = None

    eighth = -1
    cell = None
    bar = [None] * 2
    fc = [None] * 6
    fillet = [None] * 6
    raster = {'frame': -1, 'radius': 0, 'data': None}

    def __init__(self):
        self.grids = [np.zeros((self.min_size, self.min_size), dtype=bool)]
        self.grids[0][2,2] = True

    def set_bounds(self, new_bounds):
        for i in range(4):
            self.bounds[i] = range(int(new_bounds[i][0]), int(new_bounds[i][1]))

    #Return radius of snowflake #t with 10% padding
    def get_flake_radius(self, t):
        return int((len(self.grids[t]) - 4) * 1.05)

    def get_disk(self, radius, data_type):
        disk = np.zeros((radius * 2, radius * 2), dtype=data_type)

        #Generate ellipse half; optimized form of radtio^2 * (i+0.5)^2 + (j+0.5)^2 < radius^2
        for i in range(round(radius / ratio)):
            j = round((radius ** 2 - ((i + 0.5) * ratio) ** 2) ** 0.5)
            disk[radius + i, radius - j : radius + j] = 1

        #Mirror bottom half to create full disk
        disk[:radius, :] = np.flip(disk[radius:, :], 0)
        return disk

    def get_neighbor_graphic(self, neighbor, radius):
        #Highlight adjacent, colors 0 - 3
        if neighbor == 1:
            pic = np.array([[0,1,1,1,0],
                            [1,3,3,3,1],
                            [1,3,2,3,1],
                            [1,1,3,1,1],
                            [0,0,1,0,0]])

        #Highlight proximate, colors 0 - 3
        elif neighbor == 2:
            pic = np.array([[0,3,1,3,0],
                            [1,1,1,1,1],
                            [3,1,2,1,3],
                            [1,3,1,3,1],
                            [0,0,1,0,0]])

        pic_size = len(pic)

        #Create disk (filled circle) image of given radius
        disk_size = radius * 2
        disk = self.get_disk(radius, int)

        #Place disk images onto main raster according to positions given in "pic"
        raster_size = pic_size * disk_size
        raster = np.kron(pic, disk)

        for j in range(1, pic_size, 2):
            column = slice(j * disk_size, (j+1) * disk_size)
            raster[radius:-radius, column] = raster[:-disk_size, column]
            raster[:radius, column] = 0

        padding = round(raster_size * (ratio - 1) / 2)
        raster = np.pad(raster, ((0,0), (padding, padding)))
        return raster

    #Return raster of given frame
    def rasterize(self, t, radius):
        if self.raster['frame'] == t and self.raster['radius'] == radius:
            return self.raster['data']

        self.raster['frame'] = t
        self.raster['radius'] = radius
        self.raster['data'] = self.rasterize_small(t) if radius == 1 else self.rasterize_large(t, radius)

        return self.raster['data']

    #Rasterize using rectangular cells, each one a squished 2x1 image
    def rasterize_small(self, t):
        #Cell position array (A) and cell image for kronecker product
        A = self.grids[t][2:-2, 2:-2]
        cell = np.ones((2,1), dtype=bool)

        #Set up raster array
        grid_size = len(A)
        width = grid_size * 2 - 1
        raster = np.empty((width * 2, width), dtype=bool)

        #Rasterize the bottom right quarter using kronecker products, shifting every other column down 1 px
        #Two separate krons is faster than one kron followed by shifting every other column
        raster[width - 1:, grid_size - 1::2] = np.kron(A[:, ::2], cell)
        raster[width:, grid_size::2] = np.kron(A[:, 1::2], cell)[:-1, :]

        #"Unfold" raster by mirroring quarter twice
        raster[width:, :grid_size - 1] = np.flip(raster[width:, grid_size:], 1)
        raster[:width, :] = np.flip(raster[width:, :], 0)
        return raster

    #Creates pixel arrays for the bars and fillets needed to draw a smooth & rounded raster
    def gen_raster_components(self, eighth):
        self.eighth = eighth
        self.cell = self.get_disk(max(eighth * 4, 2), bool)

        #Quit now if cell is only 4x4
        if not eighth:
            return

        #Bar components between two adjacent ON cells
        self.bar[0] = np.ones((eighth*10, eighth*12), dtype=bool)
        for i in range(eighth * 4):
            self.bar[0][i, eighth * 4 + i * 2 + 1:] = False
            self.bar[0][eighth * 6 + i, :i * 2 + 1] = False
        self.bar[0][:eighth * 6, :eighth * 4] = False
        self.bar[0][eighth * 4:, eighth * 8:] = False
        self.bar[1] = np.flip(self.bar[0], 0)

        #Fillet component coordinate offsets
        self.fc[0] = [0, eighth, eighth*2, eighth*6]
        self.fc[1] = [eighth, eighth*4, eighth*6, eighth*8]
        self.fc[2] = [eighth*4, eighth*7, eighth*6, eighth*8]
        self.fc[3] = [eighth*7, eighth*8, eighth*2, eighth*6]
        self.fc[4] = [eighth*4, eighth*7, 0, eighth*2]
        self.fc[5] = [eighth, eighth*4, 0, eighth*2]

        #Fillet components for rounding off concave bends
        self.fillet[0] = np.invert(self.cell[self.fc[0][0]:self.fc[0][1], self.fc[0][2]:self.fc[0][3]])
        self.fillet[1] = np.invert(self.cell[self.fc[1][0]:self.fc[1][1], self.fc[1][2]:self.fc[1][3]])
        self.fillet[2] = np.flip(self.fillet[1], 0)
        self.fillet[3] = np.flip(self.fillet[0], 0)
        self.fillet[4] = np.flip(self.fillet[2], 1)
        self.fillet[5] = np.flip(self.fillet[1], 1)

    #Rasterize using elliptical (circular when displayed) cells, with bars and fillets between adjacent cells
    def rasterize_large(self, t, radius):
        #Generate new raster components if needed
        eighth = int(radius / 4)
        if eighth != self.eighth:
            self.gen_raster_components(eighth)

        #Get cell image dimensions
        cell_radius = max(eighth * 4, 2)
        cell_diam = cell_radius * 2

        #Load in the cell grid to be rasterized
        A = self.grids[t]
        grid_size = len(A)

        #Create empty raster. Size is double that of grid with padding removed times cell size
        raster_size = ((grid_size - 4) * 2 - 1) * cell_diam
        raster = np.empty((raster_size, raster_size), dtype=bool)

        #Generate bottom right quarter of raster using kronecker product on grid.
        half_size = int(raster_size / 2)
        offset = half_size - cell_radius
        raster[offset:, offset:] = np.kron(A[2:-2,2:-2], self.cell)

        n = [False] * 6
        #Loop through grid columns, ignoring both padding columns on the left and 3 of the 4 padding columns on the right
        for j in range(2, grid_size - 3):
            #column number of first cell pixel in raster
            cj = offset + (j - 2) * cell_diam
            parity = j % 2

            #Shift every other cell column down by half a cell
            if parity:
                raster[offset + cell_radius : -cell_radius, cj : cj + cell_diam] = raster[offset : -cell_diam, cj : cj + cell_diam]

            #Loop through grid rows, ignoring both padding rows at the top and 3 of the 4 padding rows at the bottom
            for i in range(2, grid_size - 3):
                #Row number of first cell pixel in raster
                ci = offset + (i - 2) * cell_diam + parity * cell_radius
                i_par = i + parity

                #If cell is on, copy in cell sprite and add bars between adjacent cells above & to the left
                if A[i, j]:
                    if A[i - 1, j]:
                        raster[ci - cell_radius : ci + cell_radius, cj : cj + cell_diam] = True

                    upper_left = A[i_par - 1, j - 1]
                    if upper_left and eighth:
                        raster[ci - eighth*3 : ci + eighth*7, cj - eighth*6 : cj + eighth*6] += self.bar[0]
                    elif upper_left:
                        raster[ci - 1 : ci + 3, cj - 1 : cj + 1] = True

                    lower_left = A[i_par, j - 1]
                    if lower_left and eighth:
                        raster[ci + eighth : ci + eighth*11, cj - eighth*6 : cj + eighth*6] += self.bar[1]
                    elif lower_left:
                        raster[ci + 1 : ci + 5, cj - 1: cj + 1] = True

                    continue

                #Skip fillets if cells are 4x4
                if not eighth:
                    continue

                #Put all adjacent neighbor values into a list, numbered clockwise starting from 12 o'clock
                n[0] = A[i - 1, j]
                n[1] = A[i_par - 1, j + 1]
                n[2] = A[i_par, j + 1]
                n[3] = A[i + 1, j]
                n[4] = A[i_par, j - 1]
                n[5] = A[i_par - 1, j - 1]

                #Add fillets between triplets of adjacent neighbors surrounding the 'off' cell, since they form concave bends
                for k in range(-1, 5):
                    if n[k - 1] and n[k] and n[k + 1]:
                        raster[ci + self.fc[k][0] : ci + self.fc[k][1], cj + self.fc[k][2] : cj + self.fc[k][3]] += self.fillet[k]

        #Unfold into full raster by mirroring the bottom right quarter twice
        raster[:half_size, half_size:] = np.flip(raster[half_size:, half_size:], 0)
        raster[:, :half_size] = np.flip(raster[:, half_size:], 1)

        #Invert so that 'on' cells are white
        return raster

    #Counts the near and far (adjacent and proximate) neighbors to a given cell in A
    def count_neighbors(self, A, i, j):
        #Start by accessing and combining guaranteed neighbors (independent of column number)
        near = int(A[i-1,j]) + int(A[i,j-1]) + int(A[i,j+1]) + int(A[i+1,j])
        far = int(A[i,j-2]) + int(A[i,j+2])
        above = int(A[i-1,j-1]) + int(A[i-1,j+1])
        below = int(A[i+1,j-1]) + int(A[i+1,j+1])

        #Near and far neighbors located in different rows depending on column
        if j % 2:
            near += below
            far += above + int(A[i+2,j-1]) + int(A[i+2,j+1])
        else:
            near += above
            far += below + int(A[i-2,j-1]) + int(A[i-2,j+1])

        return [near, far]

    def generate(self, t):
        #Load current generation and previous generation
        A = self.grids[t-1]
        B = self.grids[max(t-2,0)]
        size = len(A)
        prev_size = len(B)

        #Delete any future generations and get grid size
        self.grids = self.grids[:t]

        #Start with copy of current generation
        self.grids += [np.zeros((size + 2, size + 2), dtype=bool)]
        self.grids[t][2:size, 2:size] = A[2:,2:]

        #Initialize varibles used by for loops
        neigh = [0] * 4
        outermost = size - 4

        #Loop through new grid and turn on cells based on the neighbor counts in previous grids
        for i in range(2, size - 2):
            row_ncount = 0
            for j in range(2, size - 2):
                #Cell is already on, skip
                if A[i,j]:
                    continue

                #Count up near and proximate neighbors
                [neigh[0], neigh[2]] = self.count_neighbors(A, i, j)
                in_bounds = i < prev_size - 2 and j < prev_size - 2
                [neigh[1], neigh[3]] = self.count_neighbors(B, i, j) if in_bounds else [0,0]
                row_ncount += sum(neigh)

                #Check if branching rule is satisfied
                branching = self.branch and neigh[0] == 1 and (neigh[1] == 1 or neigh[2] == 0)
                
                #Alternate branching rule that creates cleaner spike intersections
                #branching = self.branch and neigh[0] == 1 and (neigh[1] == 1 or neigh[3] == 0)

                #Activate cell if on spoke or if neighbor counts are within bounds
                if branching or all([neigh[k] in self.bounds[k] for k in range(4)]):
                    self.grids[t][i,j] = True
                    outermost = max(outermost, max(i,j))

            #Finish early if no neighbors left
            if not row_ncount:
                break

        #Truncate grid to ending row & column, with 4 padding layers for next gen's neighbor counting
        flake_radius = max(i, outermost + 1)
        self.grids[t] = self.grids[t][:flake_radius + 4, :flake_radius + 4]

        #Mirror the top two rows and leftmost two columns for counting neighbors
        self.grids[t][:2, 2::2] = np.flip(self.grids[t][3:5, 2::2], 0)
        self.grids[t][:2, 3::2] = np.flip(self.grids[t][2:4, 3::2], 0)
        self.grids[t][:,:2] = np.flip(self.grids[t][:,3:5], 1)

    #Run a defined sequence of presets
    def gen_sequence(self):#, sequence)
        sequence = [[[1,2], [0,7], [0,1], [0,7], True, 13],
                    [[1,3], [0,7], [0,2], [0,7], False, 1],
                    [[1,5], [0,7], [0,7], [0,7], False, 1],
                    [[0,7], [0,7], [1,3], [1,2], False, 4]]

        total = sum([step[5] for step in sequence])
        self.__init__()
        for step in sequence:
            self.set_bounds(step[:4])
            self.spike_rule = step[4]
            for i in range(step[5]):
                frame = len(self.grids)
                self.generate(frame)
                print(f'Generating: {100 * frame/total:3.0f}%', end='\r')

class Application:
    #Delay (ms) before click-hold is repeated
    repeat_delay = 200

    #Darkest and lightest colors in 'blues' colormap
    d_blue = '#08306B'
    l_blue = '#F8FCFF'

    #Resolution radius=4 doesnt look so good
    resolutions = [2**k for k in [0, 1, 3, 4, 5]]

    #Default properties
    presets = {'Classic': [[1,2],[0,7],[0,7],[0,7], False]}
    res_num = 1
    interp = False
    flip = False
    min_radius = 2
    default_radius = 29
    zoom_radius = default_radius
    frame = 0

    #For keeping track of repeating tasks
    generating = None
    zooming = None
    advancing = None

    def __init__(self):
        #Read in presets from file and create snowflake object
        self.load_presets()
        self.sf = Snowflake()

        #Create root window
        self.root = tk.Tk()
        self.root.wm_title('Snowflake Cellular Automata')
        self.root.protocol("WM_DELETE_WINDOW", exit)

        #Estimate screen resolution using width of given font size, since winfo_fpixels sometimes gets it wrong
        font_size = 20
        font = tkfont.Font(family='Helvetica', size=font_size)
        self.font_ratio = max(font.measure("m") / font_size / 1.1, 1)

        #Create GUI input frames and pack into root
        self.create_setup()
        self.create_play()
        self.create_display_options()
        self.pack_inputs()

        if arg_split:
            self.create_flake_figure(5 * self.font_ratio) #Create flake figure of fixed size
            self.fig.show()                               #Show snowflake graphic in second window
            self.fignum = self.fig.number                 #Save figure number for checking if it gets closed
            self.btn_resize['state'] = 'disabled'         #Disable resize button since it has nothing to do
        else:
            self.root.update()                            #Update current root in order to read height in pixels for sizing snowflake frame
            self.build_GUI()                              #Builds the GUI by adding snowflake canvas (TkAgg) and re-packing input frames
            self.root.eval('tk::PlaceWindow . center')    #Re-center root window
            self.root.bind("<Key>", self.key_press)       #For binding 'R' key to the resize button

        #Start with 'Classic' preset and play buttons disabled
        self.select_preset('Classic')
        self.play_state('disabled')

    #Pack input frames into root window
    def pack_inputs(self):
        self.frm_setup.pack(side=tk.TOP, padx=10, pady=10)
        self.frm_play.pack(side=tk.BOTTOM, padx=10, pady=20)
        self.frm_display.pack(side=tk.BOTTOM, padx=10, pady=10)

    #Create snowflake's matplotlib figure with imshow plot
    def create_flake_figure(self, size):
        self.fig, self.ax = plt.subplots(figsize=(size, size), facecolor=self.d_blue)
        self.img = self.ax.imshow(self.get_image(), 'Blues', interpolation='nearest', aspect = ratio)
        self.fig.subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
        self.ax.axis('off')

    #Get window size, then clear out frames and add properly-sized snowflake frame, then re-pack input frames
    def build_GUI(self):
        #Get current window size (-10px margins) and convert into matplotlib's size units (supposedly inches)
        size = (self.root.winfo_height() - 20) / 100
        self.create_flake_figure(size)

        #Remove input frames from window to reset packing
        self.frm_setup.pack_forget()
        self.frm_play.pack_forget()
        self.frm_display.pack_forget()

        #Create oulinted frame that contains canvas for matplotlib figure
        self.frm_flake = tk.Frame(self.root, highlightbackground='grey', highlightthickness=2)
        self.canvas = FigureCanvasTkAgg(self.fig, master=self.frm_flake)
        self.canvas.get_tk_widget().pack(side=tk.TOP, fill=tk.BOTH, expand=True)

        #Pack all elements into window
        self.frm_flake.pack(side=tk.LEFT, padx=10)
        self.pack_inputs()

    #Adds the neighbor graphics, range sliders, preset menu, and generate buttons
    def create_setup(self):
        self.frm_setup = tk.LabelFrame(self.root, text='Setup', bd=5)
        self.lbl_neighbor_title = tk.Label(self.frm_setup, text='Neighbor Restrictions')

        self.neighbor_labels = ['Current Adjacent',  'Previous Adjacent', 'Current Proximate', 'Previous Proximate']
        if arg_entries:
            self.create_entries()
            self.frm_entries.grid(row=1, column=0, columnspan=3)
        else:
            self.create_sliders()
            self.cnv_neighbors.get_tk_widget().grid(row=1, column=0, columnspan=3)

        #Preset menu and branch check button
        self.frm_preset = tk.Frame(self.frm_setup)
        self.lbl_preset = tk.Label(self.frm_preset, text='Select Preset:')
        self.var_preset = tk.StringVar(self.frm_preset, 'Classic')
        self.men_preset = tk.OptionMenu(self.frm_preset, self.var_preset, *self.presets.keys(), command=self.select_preset)
        self.men_preset.config(width=22)
        self.btn_branch = tk.Checkbutton(self.frm_preset, text='Encourage branching', variable='branch', command=self.toggle_branch)

        #Generate/clear buttons
        self.btn_clear = tk.Button(master=self.frm_setup, text="Clear", font='bold', bd=4, height=2, width=9, command=self.clear)
        self.btn_gen_next = tk.Button(master=self.frm_setup, text="Generate next   ▶︎", font='bold', bd=4, height=2, width=16)
        self.btn_gen_next.bind('<ButtonPress-1>', lambda event: self.gen_next())
        self.btn_gen_next.bind('<ButtonRelease-1>', lambda event: self.stop_gen())
        self.btn_gen_all = tk.Button(master=self.frm_setup, text="Generate all   ▶︎▶︎", font='bold', bd=4, height=2, width=18, command=self.gen_all)

        #Grid placements
        self.lbl_neighbor_title.grid(row=0, column=0, columnspan=3)

        self.frm_preset.grid(row=2, column=0, columnspan=3, pady=20)
        self.lbl_preset.pack(side=tk.LEFT, padx=10)
        self.men_preset.pack(side=tk.LEFT, padx=10)
        self.btn_branch.pack(side=tk.RIGHT, padx=10)

        self.btn_clear.grid(row=3, column=0, pady=10, padx=10)
        self.btn_gen_next.grid(row=3, column=1)
        self.btn_gen_all.grid(row=3, column=2, padx=10)

    #Creates four range sliders (Matplotlib.widgets.RangeSlider) for setting neighbor bounds
    def create_sliders(self):
        #All contained in one matplotlib figure
        self.fig_neighbors = plt.figure(figsize=(4.5 * self.font_ratio , 1.5 * self.font_ratio), facecolor=self.l_blue)

        #Create the two neighbor graphics using imshow
        self.ax_neigh = [None] * 2
        for i in range(2):
            self.ax_neigh[i] = plt.axes([-0.05, 0.02 + i/2, 0.3, 0.45])
            self.ax_neigh[i].imshow(self.sf.get_neighbor_graphic(2 - i, 16), 'Blues', interpolation='bilinear', aspect = ratio)
            self.ax_neigh[i].axis('off')

        #Initialize lists for building sliders
        slider_xticks = np.arange(0, 8, 1)
        self.slider_ax = []
        self.sliders = []
        self.slider_ypositions = [0.8, 0.6, 0.3, 0.1]

        #Four total sliders
        for i in range(4):
            #Create range slider in corresponding position
            self.slider_ax += [plt.axes([0.45, self.slider_ypositions[i], 0.52, 0.1])]
            self.sliders += [RangeSlider(self.slider_ax[i], label=self.neighbor_labels[i] + '   ', valmin=0, valmax=7, valfmt='%.0f', valstep=1)]
            self.sliders[i].label.set_size(7.5 * self.font_ratio)
            self.sliders[i].valtext.set_visible(False)

            #Create slider axis. Remove major tick labels since we want the labels to be between grid marks
            self.slider_ax[i].add_artist(self.slider_ax[i].xaxis)
            self.slider_ax[i].set_xticks(np.arange(0, 8, 1))
            self.slider_ax[i].set_xticklabels('')

            #Draw integers 0 thru 6 at positions 0.5 thru 6.5 so they appear bewteen grid marks
            # e.g. when "1" is selected, the slider vals are actually [1,2]
            self.slider_ax[i].set_xticks(np.arange(0.5, 7.5, 1), minor=True)
            self.slider_ax[i].set_xticklabels(np.arange(0, 7, 1), minor=True, fontsize=7 * self.font_ratio)

            #Use grid lines rather than ticks. Grid lines are major so they align with the actual integers.
            self.slider_ax[i].grid(visible=True, which='major', color='k')
            self.slider_ax[i].tick_params(length=0, which='minor')
            self.slider_ax[i].tick_params(length=0, which='major')

        #Composite figure canvas into GUI
        self.cnv_neighbors = FigureCanvasTkAgg(self.fig_neighbors, master=self.frm_setup)

    #Creates four pairs of integer entry widgets for setting neighbor bounds
    def create_entries(self):
        self.frm_entries = tk.Frame(self.frm_setup, bg=self.l_blue)

        #Initialize lists for entries
        self.lbl_neighbor_type = [None] * 4
        self.lbl_from = [None] * 4
        self.lbl_to = [None] * 4
        padys_rel = [(6,6),(6,12),(12,6),(6,6)]
        padys = [(i[0] * self.font_ratio * 1.5, i[1] * self.font_ratio * 1.5) for i in padys_rel]

        #If initialized as [[None,None]]*4 or [[None]*2]*4, Tkinter will write var/ent objects to the same point in memory
        self.var_neighbors = [[None,None], [None,None], [None,None], [None,None]]
        self.ent_neighbors = [[None,None], [None,None], [None,None], [None,None]]

        #Four rows of entries
        for i in range(4):
            #Text labels for each row
            self.lbl_neighbor_type[i] = tk.Label(self.frm_entries, bg=self.l_blue, text=self.neighbor_labels[i] + ":")
            self.lbl_from[i] = tk.Label(self.frm_entries, bg=self.l_blue, text='from')
            self.lbl_to[i] = tk.Label(self.frm_entries, bg=self.l_blue, text='to')

            #Two entry widgets in each row
            for j in range(2):
                self.var_neighbors[i][j] = tk.StringVar(self.frm_entries, '0')
                self.ent_neighbors[i][j] = tk.Entry(self.frm_entries, width=3, textvariable=self.var_neighbors[i][j])

            #Grids placement for each row
            self.lbl_neighbor_type[i].grid(row=i, column=1, padx=(10,20), pady=padys[i])
            self.lbl_from[i].grid(row=i, column=2, padx=10, pady=padys[i])
            self.ent_neighbors[i][0].grid(row=i, column=3, padx=10, pady=padys[i])
            self.lbl_to[i].grid(row=i, column=4, padx=10, pady=padys[i])
            self.ent_neighbors[i][1].grid(row=i, column=5, padx=10, pady=padys[i])

        #If text-only GUI, return now without adding neighbor graphics
        if arg_textonly:
            return

        #Initialize neighbor graphic lists
        size = 0.7 * self.font_ratio
        self.fig_neigh = [None] * 2
        self.ax_neigh = [None] * 2
        self.cnv_neigh = [None] * 2

        #Two neighbor graphics in separate figure canvases
        for i in range(2):
            self.fig_neigh[i], self.ax_neigh[i] = plt.subplots(figsize=(size, size), facecolor=self.d_blue)
            self.ax_neigh[i].imshow(self.sf.get_neighbor_graphic(i + 1, 16), 'Blues', interpolation='bilinear', aspect = ratio)
            self.fig_neigh[i].subplots_adjust(top=1, bottom=0, right=1, left=0, hspace=0, wspace=0)
            self.ax_neigh[i].axis('off')

            self.cnv_neigh[i] = FigureCanvasTkAgg(self.fig_neigh[i], master=self.frm_entries)
            self.cnv_neigh[i].get_tk_widget().grid(row=i*2, column=0, rowspan=2, padx=10)

    #Adds the resolution, zoom, and fit buttons, along with "interpolate" and "flip" settings
    def create_display_options(self):
        self.frm_display = tk.LabelFrame(self.root, text='Display Options', bd=5)

        self.btn_reset = tk.Button(self.frm_display, text='Reset View', width=16, height=1, bd=3, command=self.reset_view)
        self.btn_resize = tk.Button(self.frm_display, text='Resize to window (R)', width=18, height=1, bd=3, command=self.rebuild_GUI)

        self.btn_resolution_minus = tk.Button(self.frm_display, text='-', width=3, height=2, bd=3, command=lambda: self.set_resolution(self.res_num - 1))
        self.lbl_resolution = tk.Label(self.frm_display, text='Resolution')
        self.btn_resolution_plus = tk.Button(self.frm_display, text='+', width=3, height=2, bd=3, command=lambda: self.set_resolution(self.res_num + 1))

        self.btn_interp = tk.Checkbutton(self.frm_display, text='Interpolate', command=self.toggle_interp)
        self.btn_flip = tk.Checkbutton(self.frm_display, text='Flip image', command=self.toggle_flip)

        self.btn_zoom_minus = tk.Button(self.frm_display, text='-', width=3, height=2, bd=3)
        self.btn_zoom_minus.bind('<ButtonPress-1>', lambda event: self.zoom(2))
        self.btn_zoom_minus.bind('<ButtonRelease-1>', lambda event: self.stop_zoom())
        self.lbl_zoom = tk.Label(self.frm_display, text='Zoom')
        self.btn_zoom_plus = tk.Button(self.frm_display, text='+', width=3, height=2, bd=3)
        self.btn_zoom_plus.bind('<ButtonPress-1>', lambda event: self.zoom(-2))
        self.btn_zoom_plus.bind('<ButtonRelease-1>', lambda event: self.stop_zoom())

        self.btn_fit_current = tk.Button(self.frm_display, text='Fit current', width=11, height=2, bd=3, command=lambda: self.set_view(self.sf.get_flake_radius(self.frame)))
        self.btn_fit_largest = tk.Button(self.frm_display, text='Fit largest', width=11, height=2, bd=3, command=lambda: self.set_view(self.sf.get_flake_radius(-1)))

        #Grid placements
        self.btn_reset.grid(row=0, column=0, pady=10, columnspan=3)
        self.btn_resize.grid(row=0, column=3, pady=10, columnspan=2)

        self.btn_resolution_minus.grid(row=1, column=0, padx=10, pady=10)
        self.lbl_resolution.grid(row=1, column=1, padx=5)
        self.btn_resolution_plus.grid(row=1, column=2, padx=10)

        self.btn_interp.grid(row=1, column=3, padx=10)
        self.btn_flip.grid(row=1, column=4, padx=(10, 20))

        self.btn_zoom_minus.grid(row=2, column=0, padx=10, pady=12)
        self.lbl_zoom.grid(row=2, column=1, padx=10)
        self.btn_zoom_plus.grid(row=2, column=2, padx=10)

        self.btn_fit_current.grid(row=2, column=3, padx=(25,10))
        self.btn_fit_largest.grid(row=2, column=4, padx=10)

    #Adds the first/prev/next/last play buttons
    def create_play(self):
        self.frm_play = tk.Frame(self.root)

        self.btn_prev = tk.Button(self.frm_play, text="◀︎", width=8, height=2, bd=3)
        self.btn_prev.bind('<ButtonPress-1>', lambda event: self.advance(-1))
        self.btn_prev.bind('<ButtonRelease-1>', lambda event: self.pause())

        self.btn_next = tk.Button(self.frm_play, text="▶︎", width=8, height=2, bd=3)
        self.btn_next.bind('<ButtonPress-1>', lambda event: self.advance(1))
        self.btn_next.bind('<ButtonRelease-1>', lambda event: self.pause())

        self.btn_first = tk.Button(self.frm_play, text="◀︎◀︎", width=6, height=2, bd=3, command=self.first)
        self.btn_last = tk.Button(self.frm_play, text="▶︎▶︎", width=6, height=2, bd=3, command=self.last)

        #Grid placements
        self.btn_first.grid(row=1, column=0, padx=10)
        self.btn_prev.grid(row=1, column=1, padx=10)
        self.btn_next.grid(row=1, column=2, padx=10)
        self.btn_last.grid(row=1, column=3, padx=10)

    #Reads in preset list from presets.txt
    def load_presets(self):
        try:
            with open('presets.txt', 'r') as pre:
                preset_list = pre.read().splitlines()
        except FileNotFoundError:
            return print('Error: ./presets.txt not found! Using "Classic" preset only')

        for entry in preset_list:
            if entry == '' or entry.startswith('#'):
                continue
            [name, data_str] = entry.rsplit(':', 1)
            data_str = data_str.strip().split(' ')
            bounds = []
            for i in range(4):
                [min, max] = data_str[i].split(',')

                bound = [0,7]
                if min.isnumeric():
                    bound[0] = int(min)
                if max.isnumeric():
                    bound[1] = int(max) + 1

                bounds += [bound]

            bounds += [bool(int(data_str[-1]))]
            self.presets.update({name.strip(): bounds})

    #Loads preset values into the sliders and updates branch flag
    def select_preset(self, selection):
        bounds = self.presets[selection]
        if arg_entries:
            for i in range(4):
                self.var_neighbors[i][0].set(bounds[i][0])
                self.var_neighbors[i][1].set(bounds[i][1] - 1)

        else:
            for i in range(4):
                #Have to do each one twice or it will occasionally not work...
                self.sliders[i].set_val(bounds[i])
                self.sliders[i].set_val(bounds[i])

            self.cnv_neighbors.draw()

        self.sf.branch = bounds[-1]
        if self.sf.branch:
            self.btn_branch.select()
        else:
            self.btn_branch.deselect()

    #Flips branch flag
    def toggle_branch(self):
        self.sf.branch = not self.sf.branch

    #Reinitializes snowflake and disables play buttons
    def clear(self):
        self.sf.__init__()
        self.zoom_radius = self.default_radius
        self.frame = 0
        self.play_state('disabled')
        self.update()

    #Sets the display state of the four play buttons
    def play_state(self, state):
        self.btn_first['state'] = state
        self.btn_prev['state'] = state
        self.btn_next['state'] = state
        self.btn_last['state'] = state

    #Stop generating upon click release
    def stop_gen(self):
        self.root.after_cancel(self.generating)
        self.generating = None

    #Read bounds from sliders or entries, updating them if needed
    def read_bounds(self):
        if not arg_entries:
            return [[int(i) for i in sld.val] for sld in self.sliders]

        bounds = [0] * 4
        for i in range(4):
            pair = [int(float(j.get())) for j in self.var_neighbors[i]]
            pair[0] = min(max(pair[0], 0), 6)
            pair[1] = min(max(pair[1], pair[0]), 6) + 1
            bounds[i] = pair

            self.var_neighbors[i][0].set(pair[0])
            self.var_neighbors[i][1].set(pair[1] - 1)

        return bounds

    #Load neighbor bounds into snowflake and generate next frame. Allows for click repeat
    def gen_next(self):
        bounds = self.read_bounds()

        self.sf.set_bounds(bounds)
        self.frame += 1
        self.sf.generate(self.frame)
        self.zoom_radius = max(self.sf.get_flake_radius(self.frame), self.zoom_radius)

        self.play_state('normal')
        self.update()

        #Initial click repeat delay of 200ms
        delay = 1
        if not self.generating:
            delay = self.repeat_delay

        self.generating = self.root.after(delay, self.gen_next)

    #Generate snowflake frames until it reaches the edge of the canvas (only displays final)
    def gen_all(self):
        bounds = self.read_bounds()
        self.sf.set_bounds(bounds)

        flake_radius = self.sf.get_flake_radius(self.frame)
        while flake_radius < self.zoom_radius:
            self.frame += 1
            self.sf.generate(self.frame)
            flake_radius = self.sf.get_flake_radius(self.frame)
            print(f'Generating: {100 * flake_radius/self.zoom_radius:3.0f}%', end='\r')

        self.frame = len(self.sf.grids) - 1
        self.play_state('normal')
        self.update()

    def reset_view(self):
        self.interp = False
        self.btn_interp.deselect()
        self.img._interpolation = 'nearest'

        self.flip = False
        self.btn_flip.deselect()

        self.set_resolution(1)
        self.set_view(self.default_radius)

    #Rebuild GUI to current window size
    def rebuild_GUI(self):
        self.set_resolution(max(self.res_num, 1))
        self.frm_flake.pack_forget()
        self.build_GUI()

    #Binding "R" to resize function
    def key_press(self, key):
        if key.char.lower() == 'r':
            self.rebuild_GUI()

    #Set resolution of snowflake given index of self.resolutions
    def set_resolution(self, new_res):
        self.res_num = new_res % len(self.resolutions)

        if self.res_num == 0:
            self.btn_resolution_minus['state'] = 'disabled'
        else:
            self.btn_resolution_minus['state'] = 'normal'

        if self.res_num == len(self.resolutions) - 1:
            self.btn_resolution_plus['state'] = 'disabled'
        else:
            self.btn_resolution_plus['state'] = 'normal'

        self.update()

    #Stops self.zoom() repitition
    def stop_zoom(self):
        self.root.after_cancel(self.zooming)
        self.zooming = None

    def zoom(self, direction):
        self.zoom_radius = max(self.min_radius, self.zoom_radius + direction)

        if self.zoom_radius == self.min_radius:
            self.btn_zoom_plus['state'] = 'disabled'
        else:
            self.btn_zoom_plus['state'] = 'normal'
        self.update()

        #Initial click repeat delay of 200ms
        delay = 1
        if not self.zooming:
            delay = self.repeat_delay

        self.zooming = self.root.after(delay, self.zoom, direction)

    #Sets zoom radius to specified value and re-enables zoom plus button if needed
    def set_view(self, size):
        self.zoom_radius = size
        self.btn_zoom_plus['state'] = 'normal'
        self.update()

    #Toggle 90-degree flip (actually transpose) for snowflake imshow
    def toggle_flip(self):
        self.flip = not self.flip
        self.update()

    #Toggle bicubic interpolation in snowflake imshow
    def toggle_interp(self):
        self.interp = not self.interp

        if self.interp:
            self.img._interpolation = 'bicubic'
        else:
            self.img._interpolation = 'nearest'
        self.update()

    #Stops self.advance() repitition
    def pause(self):
        self.root.after_cancel(self.advancing)
        self.advancing = None

    #Advances frame with repeatability
    def advance(self, direction):
        #Advance frame, cycling back to other end if needed
        self.frame = (self.frame + direction) % len(self.sf.grids)
        self.update()

        #Initial click repeat delay of 200ms, also applied when viewing final frame
        delay = 1
        if not self.advancing or self.frame == len(self.sf.grids) - 1:
            delay = self.repeat_delay

        self.advancing = self.root.after(delay, self.advance, direction)

    def first(self):
        self.frame = 0
        self.update()

    def last(self):
        self.frame = len(self.sf.grids) - 1
        self.update()

    def update(self):
        self.img.set_data(self.get_image())

        if not arg_split:
            return self.canvas.draw()

        if not plt.fignum_exists(self.fignum):
            new_fig = plt.figure(figsize=(5 * self.font_ratio, 5 * self.font_ratio))
            new_manager = new_fig.canvas.manager
            new_manager.canvas.figure = self.fig
            self.fig.set_canvas(new_manager.canvas)

        self.fig.canvas.draw()
        self.fig.show()

    def get_image(self):
        cell_radius = self.resolutions[self.res_num]
        raster = self.sf.rasterize(self.frame, cell_radius)

        cell_size = cell_radius * 2
        rast_size = len(raster)
        img_height = (self.zoom_radius * 2 + 1) * cell_size

        pad_x_div = 2 if self.res_num else 4
        pad_x = round((img_height * ratio - rast_size) / pad_x_div)
        pad_y = round((img_height - rast_size) / 2)

        if pad_y < 0:
            raster = raster[-pad_y:pad_y, :]
            pad_y = 0
        if pad_x < 0:
            raster = raster[:, -pad_x:pad_x]
            pad_x = 0

        raster = np.pad(raster, ((pad_y,pad_y), (pad_x,pad_x)))
        self.zoom_radius = int((len(raster)/cell_size - 1) / 2)
        if self.flip:
            raster = np.transpose(raster)

        #Invert so that ON is white
        return np.invert(raster)

#Create instance and run tkinter main loop
app = Application()
tk.mainloop()
