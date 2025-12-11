import Equations_finite as Ef
import Animation as An
import PlateMaterial as Pl
import In_time_animation as ItAn

plate = Pl.Plate(1, 1, 0.01, Pl.Material.Steel())

Static = Ef.Static_Euler_Bernoulli(plate, 1000, 50)
Static.Solve()

Dinamic_in_time = ItAn.In_time_Dinamic_Finite(plate, 50, 0.5, 1000)
Dinamic_in_time.Y0[:Dinamic_in_time.N_x] = Static.w
Dinamic_in_time.Create_render(5)
Dinamic_in_time.solve()

Dinamic = Ef.Euler_Bernoulli_Dinamic_Finite(plate, 50, 0.5, 50)
Dinamic.solve()
animation = An.Gif_Animation(Dinamic, 50, "giff.gif", 10)
animation.save_gif()
