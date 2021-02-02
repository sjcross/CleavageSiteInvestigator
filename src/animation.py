import math
import tkinter as tk
import time


win_w=800
win_h=600
fps = 30
step_size = 0.01

bubble_outline = "blue"
bubble_width = 3
bubble_w = 0.05
bubble_h = 0.03

dna_font = ("Ariel",12,"bold")
dna_width = 3

cass_seq_fill = "cyan"

ref_seq_fill = "magenta"

test_seq_fill = "black"

heading_font = ("Ariel",18,"bold")
heading_y = 0.05

sub_heading_font = ("Ariel", 14, "")
sub_heading_y = 0.1

cass_x1 = 0.3
cass_y1 = 0.25
cass_x2 = 0.7
cass_y2 = 0.35

ref_x1 = 0.2
ref_y1 = 0.5
ref_x2 = 0.8
ref_y2 = 0.6

test_x1 = 0.1
test_y1 = 0.75
test_x2 = 0.9
test_y2 = 0.85


# Calculating some useful values
ref_xc = ref_x1 + (ref_x2-ref_x1)/2
ref_yc = ref_y1 + (ref_y2-ref_y1)/2
cass_xc = cass_x1 + (cass_x2-cass_x1)/2
cass_yc = cass_y1 + (cass_y2-cass_y1)/2
test_xc = test_x1 + (test_x2-test_x1)/2
test_yc = test_y1 + (test_y2-test_y1)/2

# The main window of the animation
def create_window():
  window = tk.Tk()
  window.title("Animation")
  # Uses python 3.6+ string interpolation
  window.geometry(f'{win_w}x{win_h}')
  return window

  # Create a canvas for animation and add it to main window
def create_canvas(window):
  canvas = tk.Canvas(window)
  canvas.configure(bg="white")
  canvas.pack(fill="both", expand=True)
  return canvas

def create_dna(canvas, x1, y1, x2, y2, fill="black", width=1, tags=[]):
    canvas.create_line(win_w*x1, win_h*y1, win_w*x2, win_h*y1, fill=fill, width=width, tags=tags)
    canvas.create_line(win_w*x1, win_h*y2, win_w*x2, win_h*y2, fill=fill, width=width, tags=tags)

def create_text(canvas, x, y, text, font, tags=[]):
    canvas.create_text(win_w*x, win_h*y, text=text, font=font, tags=tags)

def create_bubble(canvas, x1, y1, x2, y2, bubble_outline="black", bubble_width=1, dna_fill="black", dna_width=1, tags=[]):
    yc = y1 + (y2-y1)/2
    canvas.create_line(win_w*x1, win_h*yc, win_w*x2, win_h*yc, fill=dna_fill, width=dna_width, tags=tags)
    canvas.create_rectangle(win_w*x1, win_h*y1, win_w*x2, win_h*y2,outline=bubble_outline,width=bubble_width,tags=tags)

def animate_pause(window,duration):
    n_frames = duration*fps
    window.update() # Only needs doing once, since nothing changes

    for f in range(n_frames):
        time.sleep(1/fps)
        
def animate_motion(window,canvas, x1, y1, x2, y2, step_size, tag):
    n_frames = math.ceil(math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))/step_size)
    # n_frames = duration*fps
        
    dx = win_w*(x2-x1)/n_frames
    dy = win_h*(y2-y1)/n_frames

    for f in range(n_frames):   
        canvas.move(tag,dx,dy)
        window.update()
        time.sleep(1/fps)


# Create the animaton
window = create_window()
canvas = create_canvas(window)

# Displaying DNA sequences
create_dna(canvas, ref_x1, ref_y1, ref_x2, ref_y2, fill=ref_seq_fill, width=dna_width, tags=["dna_lines","ref_seq"])
create_text(canvas, 0.5, ref_yc, "Reference sequence", dna_font, tags=["dna_text","ref_seq","text"])

create_dna(canvas, cass_x1, cass_y1, cass_x2, cass_y2, fill=cass_seq_fill, width=dna_width, tags=["dna_lines","cass_seq"])
create_text(canvas, 0.5, cass_yc, "Cassette sequence", dna_font, tags=["dna_text","cass_seq","text"])

create_dna(canvas, test_x1, test_y1, test_x2, test_y2, fill=test_seq_fill, width=dna_width, tags=["dna_lines","test_seq"])
create_text(canvas, 0.5, test_yc, "Test sequence", dna_font, tags=["dna_text","test_seq","text"])

# Displaying text
create_text(canvas, 0.5, heading_y, "Finding cassette start in test sequence", heading_font, tags=["heading","text"])
create_text(canvas, 0.5, sub_heading_y, "Extracting first N nucleotides of cassette", sub_heading_font, tags=["sub_heading"])
animate_pause(window,1)

# Showing cassette search bubble
create_bubble(canvas,cass_x1,cass_y1-bubble_h/2,cass_x1+bubble_w,cass_y1+bubble_h/2,bubble_outline=bubble_outline,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["bubble","cass_seq"])
create_bubble(canvas,cass_x1,cass_y1-bubble_h/2,cass_x1+bubble_w,cass_y1+bubble_h/2,bubble_outline=bubble_outline,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["static_bubble","cass_seq"])
animate_pause(window,1)

# Move search bubble to test sequence
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Identifying matches for cassette start in test sequence", sub_heading_font, tags=["sub_heading"])
animate_motion(window,canvas,cass_x1,cass_y1,test_x1,test_y1,step_size,"bubble")

# Scan bubble along test sequence, highlighting matching sequences
match_pos_1 = 0.3
animate_motion(window,canvas,0.12,test_y1,match_pos_1,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_1-bubble_h/2,test_y1-bubble_h/2,match_pos_1+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_fail"])
create_text(canvas, match_pos_1, test_y1-0.05, "95% match", dna_font, tags=["score_text","score_fail"])

match_pos_2 = 0.48
animate_motion(window,canvas,match_pos_1,test_y1,match_pos_2,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_2-bubble_h/2,test_y1-bubble_h/2,match_pos_2+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_fail"])
create_text(canvas, match_pos_2, test_y1-0.05, "92% match", dna_font, tags=["score_text","score_fail"])

match_pos_3 = 0.72
animate_motion(window,canvas,match_pos_2,test_y1,match_pos_3,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_3-bubble_h/2,test_y1-bubble_h/2,match_pos_3+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_pass"])
create_text(canvas, match_pos_3, test_y1-0.05, "100% match", dna_font, tags=["score_text","score_pass"])

animate_motion(window,canvas,match_pos_3,test_y1,test_x2,test_y1,step_size,"bubble")
canvas.delete("bubble")
animate_pause(window,1)

# Selecting match with highest score 
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Selecting match with highest score", sub_heading_font, tags=["sub_heading"])
canvas.delete("score_fail")
animate_pause(window,1)



name = input('Press enter to end \n')
