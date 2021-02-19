import math
import tkinter as tk
import time


win_w=800
win_h=600
fps = 30
step_size = 0.01

bubble_outline_start_sense = "orange"
bubble_outline_start_rc = "magenta"
bubble_outline_end_sense = "cyan"
bubble_outline_end_rc = "green"
bubble_width = 3
bubble_w = 0.05
bubble_h = 0.03

dna_font = ("Ariel",12,"bold")
dna_width = 3

dna_label_font = ("Ariel", 12, "")
dna_label_gap = 0.02

cass_seq_fill = "blue"
cass_x1 = 0.3
cass_y1 = 0.25
cass_x2 = 0.7
cass_y2 = 0.35

ref_seq_fill = "black"
ref_x1 = 0.2
ref_y1 = 0.5
ref_x2 = 0.8
ref_y2 = 0.6

test_seq_fill = "red"
test_x1 = 0.1
test_y1 = 0.75
test_x2 = 0.9
test_y2 = 0.85

heading_font = ("Ariel",18,"bold")
heading_y = 0.05

sub_heading_font = ("Ariel", 14, "")
sub_heading_y = 0.1

note_font = ("Ariel", 12, "")
note_y = 0.95

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

def create_dna_labels(canvas, x1, y1, x2, y2, font, gap, tags=[]):
    canvas.create_text(win_w*(x1-gap), win_h*y1, text="5'", font=font, tags=tags)
    canvas.create_text(win_w*(x2+gap), win_h*y2, text="5'", font=font, tags=tags)
    canvas.create_text(win_w*(x1-gap), win_h*y2, text="3'", font=font, tags=tags)
    canvas.create_text(win_w*(x2+gap), win_h*y1, text="3'", font=font, tags=tags)

def create_text(canvas, x, y, text, font, tags=[]):
    canvas.create_text(win_w*x, win_h*y, text=text, font=font, tags=tags)

def create_bubble(canvas, x1, y1, x2, y2, bubble_outline="black", bubble_width=1, dna_fill="black", dna_width=1, tags=[]):
    yc = y1 + (y2-y1)/2
    canvas.create_line(win_w*x1, win_h*yc, win_w*x2, win_h*yc, fill=dna_fill, width=dna_width, tags=tags)
    canvas.create_line(win_w*x1, win_h*y1, win_w*x2, win_h*y1, fill=bubble_outline, width=bubble_width, tags=tags)
    canvas.create_line(win_w*x2, win_h*y1, win_w*x2, win_h*y2, fill=bubble_outline, width=bubble_width, tags=tags)
    canvas.create_line(win_w*x2, win_h*y2, win_w*x1, win_h*y2, fill=bubble_outline, width=bubble_width, tags=tags)
    canvas.create_line(win_w*x1, win_h*y2, win_w*x1, win_h*y1, fill=bubble_outline, width=bubble_width, tags=tags)
    
def animate_pause(window,n_frames):
    window.update() # Only needs doing once, since nothing changes

    for f in range(n_frames):
        time.sleep(1/fps)
        
def animate_motion(window,canvas, x1, y1, x2, y2, step_size, tag):
    n_frames = math.ceil(math.sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1))/step_size)
            
    dx = win_w*(x2-x1)/n_frames
    dy = win_h*(y2-y1)/n_frames

    for f in range(n_frames):   
        canvas.move(tag,dx,dy)
        window.update()
        time.sleep(1/fps)

def animate_rotation(window, canvas, cx, cy, angle_degs, n_frames, tag):
    step_rads = (angle_degs * math.pi / 180)/n_frames

    ids = canvas.find_withtag(tag)

    for f in range(n_frames):   
        for id in ids:
            coords = canvas.coords(id)
            if len(coords) == 2:    
                [ox1, oy1] = canvas.coords(id)
                (rx1, ry1) = _point_rotate(ox1,oy1,win_w*cx,win_h*cy,step_rads)
                canvas.coords(id, rx1, ry1)

            elif len(coords) == 4:
                [ox1, oy1, ox2, oy2] = canvas.coords(id)
                (rx1, ry1) = _point_rotate(ox1,oy1,win_w*cx,win_h*cy,step_rads)
                (rx2, ry2) = _point_rotate(ox2,oy2,win_w*cx,win_h*cy,step_rads)
                canvas.coords(id, rx1, ry1, rx2, ry2)
        
        window.update()
        time.sleep(1/fps)

def _point_rotate(ox, oy, cx, cy, angle):
    rx = cx + math.cos(angle) * (ox - cx) - math.sin(angle) * (oy - cy)
    ry = cy + math.sin(angle) * (ox - cx) + math.cos(angle) * (oy - cy)

    return (rx, ry)


# Create the animaton
window = create_window()
canvas = create_canvas(window)
animate_pause(window,10)

# Displaying DNA sequences
create_text(canvas, 0.5, heading_y, "Loading cassette, reference and test sequences", heading_font, tags=["heading","text"])

create_dna(canvas, cass_x1, cass_y1, cass_x2, cass_y2, fill=cass_seq_fill, width=dna_width, tags=["dna_lines","cass_seq_dna"])
create_dna_labels(canvas, cass_x1, cass_y1, cass_x2, cass_y2, dna_label_font, dna_label_gap, tags=["dna_lines","cass_seq_dna"])
create_text(canvas, 0.5, cass_yc, "Cassette sequence (sense)", dna_font, tags=["dna_text","cass_seq_text","text"])

create_dna(canvas, ref_x1, ref_y1, ref_x2, ref_y2, fill=ref_seq_fill, width=dna_width, tags=["dna_lines","ref_seq_dna"])
create_dna_labels(canvas, ref_x1, ref_y1, ref_x2, ref_y2, dna_label_font, dna_label_gap, tags=["dna_lines","ref_seq_dna"])
create_text(canvas, 0.5, ref_yc, "Reference sequence (sense)", dna_font, tags=["dna_text","ref_seq_text","text"])

create_dna(canvas, test_x1, test_y1, test_x2, test_y2, fill=test_seq_fill, width=dna_width, tags=["dna_lines","test_seq_dna"])
create_dna_labels(canvas, test_x1, test_y1, test_x2, test_y2, dna_label_font, dna_label_gap, tags=["dna_lines","test_seq_dna"])
create_text(canvas, 0.5, test_yc, "Test sequence (sense)", dna_font, tags=["dna_text","test_seq_text","text"])
animate_pause(window,30)

# Showing cassette search bubble
canvas.delete("heading")
create_text(canvas, 0.5, heading_y, "Finding cassette start in test sequence", heading_font, tags=["heading","text"])
create_text(canvas, 0.5, sub_heading_y, "Extracting first N nucleotides of cassette", sub_heading_font, tags=["sub_heading"])
create_bubble(canvas,cass_x1,cass_y1-bubble_h/2,cass_x1+bubble_w,cass_y1+bubble_h/2,bubble_outline=bubble_outline_start_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["bubble"])
create_bubble(canvas,cass_x1,cass_y1-bubble_h/2,cass_x1+bubble_w,cass_y1+bubble_h/2,bubble_outline=bubble_outline_start_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["static_bubble","cass_seq_dna"])
animate_pause(window,30)

# Move search bubble to test sequence
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Identifying matches for cassette start in test sequence", sub_heading_font, tags=["sub_heading"])
animate_motion(window,canvas,cass_x1,cass_y1,test_x1,test_y1,step_size,"bubble")

# Scan bubble along test sequence, highlighting matching sequences
match_pos_start_1 = 0.3
animate_motion(window,canvas,0.12,test_y1,match_pos_start_1,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_start_1-bubble_w/2,test_y1-bubble_h/2,match_pos_start_1+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_start_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_fail"])
create_text(canvas, match_pos_start_1, test_y1-0.05, "95% match", dna_font, tags=["score_text","score_fail"])

match_pos_start_2 = 0.48
animate_motion(window,canvas,match_pos_start_1,test_y1,match_pos_start_2,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_start_2-bubble_w/2,test_y1-bubble_h/2,match_pos_start_2+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_start_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_pass"])
create_text(canvas, match_pos_start_2, test_y1-0.05, "100% match", dna_font, tags=["score_text","score_pass"])

match_pos_start_3 = 0.72
animate_motion(window,canvas,match_pos_start_2,test_y1,match_pos_start_3,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_start_3-bubble_w/2,test_y1-bubble_h/2,match_pos_start_3+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_start_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_fail"])
create_text(canvas, match_pos_start_3, test_y1-0.05, "92% match", dna_font, tags=["score_text","score_fail"])

animate_motion(window,canvas,match_pos_start_3,test_y1,test_x2,test_y1,step_size,"bubble")
canvas.delete("bubble")
animate_pause(window,30)

# Selecting match with highest score 
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Selecting match with highest score", sub_heading_font, tags=["sub_heading"])
canvas.delete("score_fail")
animate_pause(window,30)

# Get reverse complement of cassette
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Extracting first N nucleotides of cassette reverse complement", sub_heading_font, tags=["sub_heading"])
animate_rotation(window,canvas,0.5,cass_yc,180,30,"cass_seq_dna")
canvas.itemconfigure("cass_seq_text",text="Cassette sequence (RC)")
create_bubble(canvas,cass_x1,cass_y1-bubble_h/2,cass_x1+bubble_w,cass_y1+bubble_h/2,bubble_outline=bubble_outline_start_rc,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["bubble"])
create_bubble(canvas,cass_x1,cass_y1-bubble_h/2,cass_x1+bubble_w,cass_y1+bubble_h/2,bubble_outline=bubble_outline_start_rc,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["static_bubble","cass_seq_dna"])
animate_pause(window,30)

# Move search bubble to test sequence
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Identifying matches for RC cassette start in test sequence", sub_heading_font, tags=["sub_heading"])
animate_motion(window,canvas,cass_x1,cass_y1,test_x1,test_y1,step_size,"bubble")

# Scan bubble along test sequence, highlighting matching sequences
match_pos_start_rc_1 = 0.2
animate_motion(window,canvas,0.12,test_y1,match_pos_start_rc_1,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_start_rc_1-bubble_w/2,test_y1-bubble_h/2,match_pos_start_rc_1+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_start_rc,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_pass_temp"])
create_text(canvas, match_pos_start_rc_1, test_y1-0.05, "98% match", dna_font, tags=["score_text","score_pass_temp"])

match_pos_start_rc_2 = 0.67
animate_motion(window,canvas,match_pos_start_rc_1,test_y1,match_pos_start_rc_2,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_start_rc_2-bubble_w/2,test_y1-bubble_h/2,match_pos_start_rc_2+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_start_rc,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_fail"])
create_text(canvas, match_pos_start_rc_2, test_y1-0.05, "94% match", dna_font, tags=["score_text","score_fail"])

animate_motion(window,canvas,match_pos_start_rc_2,test_y1,test_x2,test_y1,step_size,"bubble")
canvas.delete("bubble")
animate_pause(window,30)

# Selecting match with highest score 
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Selecting RC match with highest score", sub_heading_font, tags=["sub_heading"])
canvas.delete("score_fail")
animate_pause(window,30)

# Selecting overall match with highest score 
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Selecting overall match with highest score (sense and RC)", sub_heading_font, tags=["sub_heading"])
canvas.delete("score_pass_temp")
animate_pause(window,30)
canvas.delete("score_text")


# Finding cassette end in same orientation as highest-match cassette start
canvas.delete("heading")
canvas.delete("sub_heading")
create_text(canvas, 0.5, heading_y, "Finding cassette end in test sequence", heading_font, tags=["heading","text"])
create_text(canvas, 0.5, sub_heading_y, "Extracting cassette end in same orientation as highest-match cassette start*", sub_heading_font, tags=["sub_heading"])
create_text(canvas, 0.5, note_y, "*In this case, the sense orientation", note_font, tags=["note"])
animate_rotation(window,canvas,0.5,cass_yc,180,30,"cass_seq_dna")
canvas.itemconfigure("cass_seq_text",text="Cassette sequence (sense)")

# Get last N nucleotides of RC cassette
create_bubble(canvas,cass_x2-bubble_w,cass_y1-bubble_h/2,cass_x2,cass_y1+bubble_h/2,bubble_outline=bubble_outline_end_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["bubble"])
create_bubble(canvas,cass_x2-bubble_w,cass_y1-bubble_h/2,cass_x2,cass_y1+bubble_h/2,bubble_outline=bubble_outline_end_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["static_bubble","cass_seq_dna"])
animate_pause(window,30)

# Move search bubble to test sequence
canvas.delete("sub_heading")
canvas.delete("note")
create_text(canvas, 0.5, sub_heading_y, "Identifying matches for cassette end in test sequence", sub_heading_font, tags=["sub_heading"])
animate_motion(window,canvas,cass_x2-bubble_w,cass_y1,test_x1,test_y1,step_size,"bubble")

# Scan bubble along test sequence, highlighting matching sequences
match_pos_end_1 = 0.32
animate_motion(window,canvas,0.12,test_y1,match_pos_end_1,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_end_1-bubble_w/2,test_y1-bubble_h/2,match_pos_end_1+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_end_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_fail"])
create_text(canvas, match_pos_end_1, test_y1-0.05, "95% match", dna_font, tags=["score_text","score_fail"])

match_pos_end_2 = 0.62
animate_motion(window,canvas,match_pos_end_1,test_y1,match_pos_end_2,test_y1,step_size,"bubble")
create_bubble(canvas,match_pos_end_2-bubble_w/2,test_y1-bubble_h/2,match_pos_end_2+bubble_w/2,test_y1+bubble_h/2,bubble_outline=bubble_outline_end_sense,bubble_width=bubble_width,dna_fill=cass_seq_fill,dna_width=dna_width,tags=["candidate","score_pass"])
create_text(canvas, match_pos_end_2, test_y1-0.05, "100% match", dna_font, tags=["score_text","score_pass"])

animate_motion(window,canvas,match_pos_end_2,test_y1,test_x2,test_y1,step_size,"bubble")
canvas.delete("bubble")
animate_pause(window,30)

# Selecting match with highest score 
canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Selecting match with highest score", sub_heading_font, tags=["sub_heading"])
canvas.delete("score_fail")
animate_pause(window,30)

canvas.delete("sub_heading")
create_text(canvas, 0.5, sub_heading_y, "Cassette found in test sequence", sub_heading_font, tags=["sub_heading"])
canvas.delete("score_text")
create_dna(canvas, match_pos_start_2-bubble_w/2, test_y1, match_pos_end_2+bubble_w/2, test_y2, fill=cass_seq_fill, width=dna_width, tags=["dna_lines","cass_seq_in_test_dna"])
animate_pause(window,30)
canvas.delete("bubble")

name = input('Press enter to end \n')
