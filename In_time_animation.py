import matplotlib.pyplot
import scipy.integrate
import numpy
import PlateMaterial
import time


class In_time_render():
    def __init__(self, L, dt, x, y):
        matplotlib.pyplot.ion()
        self.x = x
        self.dt = dt
        self.figure, self.axis = matplotlib.pyplot.subplots()

        self.line, = self.axis.plot(self.x, y, 'b-', linewidth=2)
        self.axis.set_xlim(0, L)

        self.axis.set_ylim(-0.3, 0.3)

        self.axis.set_xlabel('Длина пластины, м')
        self.axis.set_ylabel('Прогиб, мм')
        self.axis.grid(True, alpha=0.3)

        matplotlib.pyplot.show(block=False)
        matplotlib.pyplot.pause(0.5)

    def Update(self, Y, time):
        current_time = time
        self.axis.set_title(f'Время: {current_time:.3f} с')
        self.line.set_ydata(Y)
        self.figure.canvas.draw_idle()
        self.figure.canvas.flush_events()


class In_time_Dinamic_Finite():
    """
    Решает уравнение колебаний Эйлера - бернулли, методом конечных элементов
    Создаёт симуляцию в реальном времени
    Принимает объект класса пластины, количество элементов n_x,
    Время симуляции T, Количество шагов симуляции n_t
    """

    def __init__(self, plate, n_x, T, n_t):
        self.L = plate.Get_l()
        self.N_t = n_t
        self.N_x = n_x
        self.E = plate.Get_E()
        self.I = plate.Get_I()
        self.nu = plate.Get_nu()
        self.T = T

        self.dx = self.L / (self.N_x - 1)
        self.t = numpy.linspace(0, self.T, self.N_t)

        # Первое значение - кривизна w второе производная кривизны dwdt
        self.Y0 = numpy.zeros((2 * self.N_x))
        self.result = None

    def Create_render(self, update_frequency):
        self.update_counter = 0
        self.update_frequency = update_frequency
        self.in_time_render = In_time_render(
            self.L, self.T / self.N_t, numpy.linspace(0, self.L, self.N_x), self.Y0[:self.N_x])

    def Euler_Equation(self, Y, t):
        """
        Уравнение используемое функци odeit для определения дифференциального уравнения
        Функция принимает Y - двумерный массив
        Y[0] - массив кривизны в точке, длинной N_x
        Y[1] - массив первой производной кривизны в точке, длинной N_x
        t - свободный аргумент, время по которому проводится решение
        N_x - количество точек пластины, по которым проводится решение
        dx - расстояние между точками пластины
        E - модуль юнга
        I - момент инерции пластины в продольном направлении
        nu - масса на еденицу длинны, линейная плотность пластины
        """
        dx = self.dx
        E = self.E
        I = self.I
        nu = self.nu

        w = Y[:self.N_x]
        dwdt = Y[self.N_x:]
        result = numpy.zeros((2 * self.N_x))

        self.update_counter += 1
        if self.update_counter % self.update_frequency == 0:
            self.in_time_render.Update(w, t)
            time.sleep(0.001)

        dwdt[0] = 0
        dwdt[self.N_x - 1] = 0
        result[:self.N_x] = dwdt  # cкорость точек

        dw2dt2 = numpy.zeros((self.N_x))  # ускорение точек

        for i in range(2, self.N_x - 2):

            # пересчёт четвёртой производной приближенным методом
            d4w_dx4 = (w[i - 2] - 4*w[i - 1] + 6*w[i] -
                       4*w[i + 1] + w[i + 2]) / dx**4

            dw2dt2[i] = -E * I / nu * d4w_dx4 / 1000  # в мм / с^2

        # Обработка крайних точек
        dw2dt2[self.N_x - 1] = -E * I / nu * \
            (w[self.N_x - 3] - 2*w[self.N_x - 2] +
             w[self.N_x - 1]) / dx**4 / 1000

        dw2dt2[self.N_x - 2] = -E * I / nu * \
            (w[self.N_x - 4] - 4*w[self.N_x - 3] + 5 *
             w[self.N_x - 2] - 2*w[self.N_x - 1]) / dx**4 / 1000

        dw2dt2[1] = -E * I / nu * \
            (w[3] - 4*w[2] + 5 *
             w[1] - 2*w[0]) / dx**4 / 1000

        dw2dt2[0] = 0

        result[self.N_x:] = dw2dt2

        return result

    def solve(self):
        """
        Решает дифференциальное уравнение, возвращает массивы значания кривизны в точках.
        """
        print("Started solving dinamic")
        result = scipy.integrate.odeint(
            self.Euler_Equation, self.Y0, self.t, mxstep=100000, rtol=1e-8, atol=1e-8)

        print("Ended solving dinamic")
        self.result = result
