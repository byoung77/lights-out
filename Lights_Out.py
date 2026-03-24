import argparse
import tkinter as tk
import tkinter.font as tkfont
from tkinter import ttk
from tkinter import messagebox
from PIL import Image, ImageTk
import sys
import numpy as np
import os
from itertools import product

script_dir = os.path.dirname(os.path.abspath(__file__))

class LightsOutException(Exception):
    pass

#colors dictionary
colors_dict = {2:["gray41","gold"], 
                3:["gray41","red", "gold"], 
                5:["gray41","red","green3","blue","gold"],
                7:["gray41","red","orange","green3","DeepSkyBlue2","blue violet","gold"]}

def press_pluralizer(press):
    if press > 1:
        return "presses"
    return "press"

def create_adj_plus(adj_matrix, p):
    """
    Compute a reduced form of the adjacency matrix over Z_p together with
    a solver matrix A^+ obtained from the same row operations.

    For any b in the image of A, A^+ b gives a solution to the reduced system,
    and in this application yields a valid press vector solving A x = b mod p.
    """
    dim = adj_matrix.shape[0]
    right_matrix = np.identity(dim, dtype=int)
    full_matrix = np.hstack([adj_matrix, right_matrix])
    
    pivot_row = 0
    for col in range(dim):
        pivot = None
        # find a pivot
        for r in range(pivot_row, dim):
            if full_matrix[r,col] % p != 0:
                pivot = r
                break
        if pivot is None:
            continue
    
        #swap rows if needed
        if pivot != pivot_row:
            full_matrix[[pivot_row, pivot]] = full_matrix[[pivot, pivot_row]]
            
        #scale pivot row so pivot entry becomes 1 mod p
        pivot_val = full_matrix[pivot_row, col] % p
        pivot_inv = pow(int(pivot_val), -1, p)
        full_matrix[pivot_row] = (pivot_inv * full_matrix[pivot_row]) % p

        #eliminate this column from all other rows
        for r in range(dim):
            if r != pivot_row and full_matrix[r, col] % p != 0:
                factor = full_matrix[r, col] % p
                full_matrix[r] = (full_matrix[r] - factor * full_matrix[pivot_row]) % p
                
        pivot_row += 1
        if pivot_row == dim:
            break
    
    return full_matrix[:,:dim], full_matrix[:,dim:]



#function to create button press neighbor function
def button_press_gen(i,j,rows,cols,top):
    all_nbrs = [(i,j-1), (i-1,j), (i+1,j), (i,j+1)]
    
    if top == 'simple':
        # adjust neighbors
        ret_nbrs = []
        for nbr in all_nbrs:
            if nbr[0]==-1 or nbr[0]==rows:
                continue
            elif nbr[1]==-1 or nbr[1]==cols:
                continue
            else:
                ret_nbrs.append(nbr)
                
        ret_nbrs = list(set(ret_nbrs))
        
        button_press_nbrs = lambda : ret_nbrs
        return button_press_nbrs
    
    if top == 'cylinder':
        # adjust neighbors
        ret_nbrs = []
        for nbr in all_nbrs:
            if nbr[0]==-1 or nbr[0]==rows:
                continue
            elif nbr[1]==-1:
                ret_nbrs.append((nbr[0],cols-1))
            elif nbr[1]==cols:
                ret_nbrs.append((nbr[0],0))    
            else:
                ret_nbrs.append(nbr)
        ret_nbrs = list(set(ret_nbrs))
        
        button_press_nbrs = lambda : ret_nbrs
        return button_press_nbrs
    
    if top == 'mobius':
        # adjust neighbors
        ret_nbrs = []
        for nbr in all_nbrs:
            if nbr[0]==-1 or nbr[0]==rows:
                continue
            elif nbr[1]==-1:
                ret_nbrs.append((rows-1-nbr[0],cols-1))
            elif nbr[1]==cols:
                ret_nbrs.append((rows-1-nbr[0],0))    
            else:
                ret_nbrs.append(nbr)
                
        ret_nbrs = list(set(ret_nbrs))
        
        button_press_nbrs = lambda : ret_nbrs
        return button_press_nbrs
    
    if top == 'torus':
        # adjust neighbors
        ret_nbrs = []
        for nbr in all_nbrs:
            if nbr[0]==-1:
                ret_nbrs.append((rows-1,nbr[1]))
            elif nbr[0]==rows:
                ret_nbrs.append((0,nbr[1]))
            elif nbr[1]==-1:
                ret_nbrs.append((nbr[0],cols-1))
            elif nbr[1]==cols:
                ret_nbrs.append((nbr[0],0))    
            else:
                ret_nbrs.append(nbr)
        
        ret_nbrs = list(set(ret_nbrs))
        
        button_press_nbrs = lambda : ret_nbrs
        return button_press_nbrs
    
    if top == 'klein':
        # adjust neighbors
        ret_nbrs = []
        for nbr in all_nbrs:
            if nbr[0]==-1:
                ret_nbrs.append((rows-1,nbr[1]))
            elif nbr[0]==rows:
                ret_nbrs.append((0,nbr[1]))
            elif nbr[1]==-1:
                ret_nbrs.append((rows-1-nbr[0],cols-1))
            elif nbr[1]==cols:
                ret_nbrs.append((rows-1-nbr[0],0))    
            else:
                ret_nbrs.append(nbr)
        
        ret_nbrs = list(set(ret_nbrs))
        
        button_press_nbrs = lambda : ret_nbrs
        return button_press_nbrs
    
    if top == 'projective':
        # adjust neighbors
        ret_nbrs = []
        for nbr in all_nbrs:
            if nbr[0]==-1:
                ret_nbrs.append((rows-1,cols-1-nbr[1]))
            elif nbr[0]==rows:
                ret_nbrs.append((0,cols-1-nbr[1]))
            elif nbr[1]==-1:
                ret_nbrs.append((rows-1-nbr[0],cols-1))
            elif nbr[1]==cols:
                ret_nbrs.append((rows-1-nbr[0],0))    
            else:
                ret_nbrs.append(nbr)
        
        ret_nbrs = list(set(ret_nbrs))
        
        button_press_nbrs = lambda : ret_nbrs
        return button_press_nbrs

#button class
class Button:
    def __init__(self, i, j, rows, cols, top, tot_states):
        self.row = i
        self.col = j
        self.topology = top
        self.tot_states = tot_states
        self.state = 0
        self.color_list = colors_dict[tot_states]
        self.color = self.color_list[0]
        self.call_nbrs = button_press_gen(i,j,rows,cols,top)
        
    def __call__(self):
        self.state = (self.state + 1) % self.tot_states
        self.color = self.color_list[self.state]
        return self.call_nbrs()
        
    def update_from_nbr(self):
        self.state = (self.state + 1) % self.tot_states
        self.color = self.color_list[self.state]
        
class LightsOutGrid:
    def __init__(self, rows, cols, top, tot_states):
        self.rows = rows
        self.cols = cols
        self.topology = top
        self.tot_states = tot_states
        self.init_state = None
        
        self.game_grid = []
        for i in range(rows):
            new_row = []
            for j in range(cols):
                new_row.append(Button(i,j,rows,cols,top,tot_states)) 
            self.game_grid.append(new_row)
        
        self.adj_matrix = np.identity(rows*cols, dtype=int)
        for i,j in product(range(rows), range(cols)):
            col = i*cols + j
            nbrs = self.game_grid[i][j].call_nbrs()
            for r,c in nbrs:
                row = r*cols + c
                self.adj_matrix[row, col] += 1
        self.adj_matrix %= tot_states
        
        self.adj_matrix_red, self.adj_matrix_plus = create_adj_plus(self.adj_matrix, tot_states)
        self.adj_nullity = np.sum(np.all(self.adj_matrix_red==0,axis=1))

    
    def __call__(self, i, j):
        nbrs = self.game_grid[i][j]()
        for nbr in nbrs:
            self.game_grid[nbr[0]][nbr[1]].update_from_nbr()
    
    def grid_loc_to_vec_idx(self,i,j):
            return i*self.cols + j 
        
    def vec_idx_to_grid_loc(self, idx):
        j = idx % self.cols
        i = (idx - j)//self.cols
        return i,j

    def state_vector(self):
        return np.array([btn.state for row in self.game_grid for btn in row], dtype=int)

    def is_solved(self):
        return all(state==0 for state in self.state_vector())
    
    def initialize_state(self):
        while True:
            win_vec = np.random.choice(range(self.tot_states), size=self.rows*self.cols)
            init_vec = self.adj_matrix @ win_vec
            init_vec %= self.tot_states
            if np.any(init_vec != 0):
                break
        self.init_state = init_vec.copy()

        for idx in range(self.rows*self.cols):
            s = init_vec[idx]
            i,j = self.vec_idx_to_grid_loc(idx)
            self.game_grid[i][j].state = s
            self.game_grid[i][j].color = self.game_grid[i][j].color_list[s]
            
    def reset_grid(self):
        for idx in range(self.rows*self.cols):
            s = self.init_state[idx]
            i,j = self.vec_idx_to_grid_loc(idx)
            self.game_grid[i][j].state = s
            self.game_grid[i][j].color = self.game_grid[i][j].color_list[s]
        
    def generate_soln(self):
        #get solution vector by solving A@soln_vec = -1*current_state (mod tot_states)
        neg_state_vec = (-1*self.state_vector()) % self.tot_states
        soln_vec = (self.adj_matrix_plus @ neg_state_vec) % self.tot_states
        
        #assemble hint string
        if self.adj_nullity==0:
            hint_str = "The Solution:\n======================"
        else:
            hint_str = "One Possible Solution:\n======================"
        for idx in range(len(soln_vec)):
            num_press = soln_vec[idx]
            if num_press != 0:
                i, j = self.vec_idx_to_grid_loc(idx)
                hint_str += f"\nButton in row={i+1:>2}, column={j+1:>2}: {num_press} {press_pluralizer(num_press)}"
        return hint_str
    
def start_lights_out_game(root, rows, cols, s, top):
    """
    # Old Argument Parsing for CLI version
    parser = argparse.ArgumentParser(description="Simulates the Lights Out game!  Default game is a 5X5 grid with 2 states.")
    parser.add_argument("-r", dest="rows", default=5, type=int, help="sets the number of rows (between 2 and 15)")
    parser.add_argument("-c", dest="cols", default=5, type=int, help="sets the number of columns (between 2 and 15)")
    parser.add_argument("-s", dest="states", default=2, type=int, help="sets the number of states (one of 2, 3, 5, or 7)")
    parser.add_argument("-t", dest="topology", default='simple', type=str, help="sets the board topology (types: simple (default), cylinder, mobius, torus, klein, or projective)")
    args = parser.parse_args()
    rows = args.rows
    cols = args.cols
    s = args.states
    top = args.topology
    
    #check passed arguments
    if rows > 15:
        raise LightsOutException("Given row size is too large (max is 15).")

    if rows < 2:
        raise LightsOutException("Given row size is too small (min is 2).")

    if cols > 15:
        raise LightsOutException("Given column size is too large (max is 15).")

    if cols < 2:
        raise LightsOutException("Given column size is too small (min is 2).")
        
    if s not in [2,3,5,7]:
        raise LightsOutException("Total number of states must be 2, 3, 5, or 7.")
    
    if top not in ["simple", "cylinder", "mobius", "torus", "klein", "projective"]:
        raise LightsOutException("Unrecognized topology (options are simple, cylinder, mobius, torus, klein, or projective).")
    
    """
    
    for widget in root.winfo_children():
        widget.destroy()

     
    game = LightsOutGrid(rows,cols,top,s)
    game.initialize_state()
    
    btn_height = 2
    btn_width = 3
    
    if rows <= 10:
        btn_height = 3
        btn_width = 5
    
    button_array = []
    for i in range(rows):
        new_row = []
        for j in range(cols):
            new_button = tk.Button(root, bg=game.game_grid[i][j].color, command = lambda i1 = i, j1 = j: update_grid(i1, j1))
            new_button.config(height = btn_height, width = btn_width)
            new_button.config(activebackground=new_button.cget("background"))
            new_button.grid(row = i, column = j, padx = 3, pady = 3)
            new_row.append(new_button)
        button_array.append(new_row)
        
    def mouseClick(event):
        game_launcher(root)
        
    def close_game():
        root.destroy()
    
    def menu_action():
        game_launcher(root)
    
    def redraw_grid():
        if game.is_solved():
            for widget in root.winfo_children():
                widget.destroy()
                
            root.geometry("800x600")
            root.win_icon = tk.PhotoImage(file=os.path.join(script_dir, "images/lightsout.png"))
            root.iconphoto(True, root.win_icon)
            load = Image.open(os.path.join(script_dir, "images/youwin.jpg"))
            render = ImageTk.PhotoImage(load)
            img = tk.Label(root, image=render)
            img.image = render
            img.bind( "<Button>", mouseClick ) 
            img.place(x=0, y=0)
        
        else:
            for i,j in product(range(rows),range(cols)):
                button_array[i][j].config(bg=game.game_grid[i][j].color)
                button_array[i][j].config(activebackground=button_array[i][j].cget("background"))
    
    # reset function
    def reset():
        game.reset_grid()
        redraw_grid()
    
    # new game function
    def new_game():
        game.initialize_state()
        redraw_grid()
        
    def update_grid(i,j):
        game(i,j)
        redraw_grid()
        
    def show_hint():
        hint_win = tk.Toplevel(root)
        hint_win.title("Hint")
        hint_win.resizable(True, True)

        # Frame to hold text + scrollbar
        frame = tk.Frame(hint_win)
        frame.pack(fill="both", expand=True, padx=10, pady=10)

        # Text widget
        txt = tk.Text(
            frame,
            wrap="none",
            font=("Courier New", 10)
        )
        txt.pack(side="left", fill="both", expand=True)

        # Vertical scrollbar
        scrollbar = tk.Scrollbar(frame, orient="vertical", command=txt.yview)
        scrollbar.pack(side="right", fill="y")

        txt.config(yscrollcommand=scrollbar.set)

        # Insert content
        txt.insert("1.0", game.generate_soln())
        txt.config(state="disabled")

        # Close button
        close_hint = tk.Button(hint_win, text="Close", command=hint_win.destroy)
        close_hint.pack(pady=(0,10))
    #reset and new buttons
    resetButton = tk.Button(root, text = "Reset", command = reset)
    resetButton.config(width = btn_width)
    resetButton.grid(row = rows+1, column = 0, padx = 3, pady = 3)
    newButton = tk.Button(root, text = "New", command = new_game)
    newButton.config(width = btn_width)
    newButton.grid(row = rows+1, column = 1, padx = 3, pady = 3)
    
    hint_col = cols-2
    if hint_col < 2:
        hint_col = 2
        
    menu_col = cols-1
    if menu_col < 3:
        menu_col = 3
        
    hintButton = tk.Button(root, text = "Hint", command = show_hint)
    hintButton.config(width = btn_width)
    hintButton.grid(row = rows+1, column = hint_col, padx = 3, pady = 3)
    menuButton = tk.Button(root, text = "Menu", command = menu_action)
    menuButton.config(width =btn_width)
    menuButton.grid(row = rows+1, column = menu_col, padx = 3, pady = 3)
    
    root.update_idletasks()

    w = root.winfo_reqwidth()
    h = root.winfo_reqheight()

    screen_w = root.winfo_screenwidth()
    screen_h = root.winfo_screenheight()

    x = max((screen_w - w) // 2, 0)
    y = max((screen_h - h) // 2, 0)

    root.geometry(f"{w}x{h}+{x}+{y}")
    

    
def game_launcher(root):
    for widget in root.winfo_children():
        widget.destroy()
        
    label_font = tkfont.Font(size=12)
    widget_font = tkfont.Font(size=12)
    button_font = tkfont.Font(size=12, weight="bold")
    
    def update_topology_image(event=None):
        photo = loaded_images[top_var.get()]
        image_label.config(image=photo, text="")
        image_label.image = photo
        
    def close_action():
        root.destroy()
        
    def play():
        #get data from entry window
        rows = rows_var.get()
        cols = cols_var.get()
        states = states_var.get()
        top = top_var.get()
        
        #launch game
        start_lights_out_game(root, rows, cols, states, top)
    
    #row and column choices  
    rows_var = tk.IntVar(value=5)
    cols_var = tk.IntVar(value=5)
    states_var = tk.IntVar(value=2)
    top_var = tk.StringVar(value='simple')

    size_choices = [str(i) for i in range(2, 16)]
    state_choices = [2,3,5,7]
    top_choices = ["simple", "cylinder", "mobius", "torus", "klein", "projective"]
    top_images = {"simple":os.path.join(script_dir, "images/simple_grid.png"), "cylinder":os.path.join(script_dir, "images/cylinder_grid.png"), 
					"mobius":os.path.join(script_dir, "images/mobius_grid.png"), "torus":os.path.join(script_dir, "images/torus_grid.png"),
					"klein":os.path.join(script_dir, "images/klein_grid.png"), "projective":os.path.join(script_dir, "images/projective_grid.png")}
    loaded_images = {k: ImageTk.PhotoImage(Image.open(v)) for k, v in top_images.items()}

    tk.Label(root, text="Rows", font=label_font).grid(row=0, column=0, padx=5, pady=5)
    rows_menu = ttk.Spinbox(root, textvariable=rows_var, values=size_choices, state="readonly", width=4,  font=widget_font)
    rows_menu.grid(row=0, column=1, padx=5, pady=5)

    tk.Label(root, text="Columns", font=label_font).grid(row=1, column=0, padx=5, pady=5)
    cols_menu = ttk.Spinbox(root, textvariable=cols_var, values=size_choices, state="readonly", width=4,  font=widget_font)
    cols_menu.grid(row=1, column=1, padx=5, pady=5)
    
    tk.Label(root, text="States", font=label_font).grid(row=0, column=3, padx=5, pady=5)
    state_menu = ttk.Combobox(root, textvariable=states_var, values=state_choices, state="readonly", width=4,  font=widget_font)
    state_menu.grid(row=0, column=4, padx=5, pady=5)
    
    tk.Label(root, text="Topology", font=label_font).grid(row=1, column=3, padx=5, pady=5)
    top_menu = ttk.Combobox(root, textvariable=top_var, values=top_choices, state="readonly", width=9,  font=widget_font)
    top_menu.grid(row=1, column=4, padx=5, pady=5)
    
    top_menu.current(0)
    
    image_label = tk.Label(root)
    image_label.grid(row=2, column=0, columnspan=5, padx=5, pady=10)

    top_menu.bind("<<ComboboxSelected>>", update_topology_image)
            
    playButton = tk.Button(root, text = "Play", font=button_font, command = play)
    playButton.config(width = 5)
    playButton.grid(row = 3, column = 0, padx = 5, pady = 5)
    
    closeButton = tk.Button(root, text = "Exit", font=button_font, command = close_action)
    closeButton.config(width = 5)
    closeButton.grid(row = 3, column = 4, padx = 5, pady = 5)
    
    root.update_idletasks()
    update_topology_image()
    root.geometry("") 
    root.update()
    
    
    

if __name__ == "__main__":
    np.random.seed(os.getpid())
    
    root = tk.Tk(className="LightsOut")
    root.wm_title("Lights Out")
    
    icon_path = os.path.join(script_dir, "images/lightsout.png")
    root.icon_img = tk.PhotoImage(file=icon_path)
    root.iconphoto(True, root.icon_img)
    
    game_launcher(root)
    
    #Main Loop
    root.mainloop()
    

                

        
    

