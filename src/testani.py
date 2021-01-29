import tkinter as yk
import tkinter as tk

window = tk.Tk()
window.title("Animation")

canvas = tk.Canvas(window)
canvas.configure(bg="white")
canvas.pack(fill="both", expand=True)

canvas.create_line(40,40,40,90,tags=["dna","left"])
canvas.create_line(60,40,60,90,tags=["dna","right"])
canvas.create_line(40,100,40,190,tags=["dna","left"])
canvas.create_line(60,100,60,190,tags=["dna","right"])

name = input('Press enter to move DNA \n')
canvas.move("dna",50,0)
name = input('Press enter to move right DNA only \n')
canvas.move("right",50,0)
name = input('Press enter to end \n')
