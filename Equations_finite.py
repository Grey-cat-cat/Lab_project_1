import matplotlib.pyplot
import scipy.integrate
import numpy


class Static_plate():
    '''
    Создаёт прогиб пластины, при постоянной силе, воздействующей на центр
    Не поддерживает анимацию, но может вывеестиграфик
    Принимает на вход объект plate класаа Plate, количество точек n_x и силу F
    '''

    def __init__(self, plate, F, n_x):
        self.L = plate.Get_l()
        self.n_x = n_x
        self.E = plate.Get_E()
        self.I = plate.Get_I()
        self.F = F
        self.dx = self.L / (self.n_x - 1)

        self.x = numpy.linspace(0, self.L, self.n_x)

        self.Y0 = numpy.zeros((2))
        # 2 мерный массив, содержит значения кривизны и её производных в каждой точке
        self.result = None
        # Значение кривизны в каждой точке, одномерный массив
        self.w = numpy.zeros((self.n_x))

    def calculete(self, Y, x):
        """
        Решение дифференциального уравнения для определения прогиба пластины, под действием постоянной силы
        """
        w, u1 = Y

        return [u1, -(-x) * self.F / (2 * self.E * self.I) * (1) ** (3/2)]

    def Solve(self):
        '''
        Записывает результат дифференциального уравнения в ячейку объекта self.result
        Записывает значения прогиба в одномерный массив self.w
        '''
        self.result = scipy.integrate.odeint(
            self.calculete, self.Y0, self.x, mxstep=100000, rtol=1e-8, atol=1e-10)

        for i in range(self.n_x):
            self.w[i] = self.result[i][0]

    def Get_render(self):
        """
        Выводит график значений координаты от высоты для каждой точки пластины
        """
        array_y = numpy.zeros((2 * self.n_x))
        array_x = numpy.zeros((2 * self.n_x))

        array_x = numpy.zeros((2 * self.n_x))
        for i in range(self.n_x):  # Пересчёт координаты x
            array_x[self.n_x + i] = self.dx * i
            array_x[self.n_x - 1 - i] = self.dx * i
            array_y[self.n_x + i] = self.w[i]
            array_y[self.n_x - 1 - i] = self.w[i]

        matplotlib.pyplot.plot(
            array_x, array_y, 'b-', linewidth=2)
        matplotlib.pyplot.show()

    def hight(self):
        """Возвращает максимальныйй прогиб"""
        return self.w[self.n_x - 1]


class Static_Euler_Bernoulli():

    def __init__(self, plate, F, n_x):
        self.L = plate.Get_l()
        self.n_x = n_x
        self.E = plate.Get_E()
        self.I = plate.Get_I()
        self.q = F / plate.Get_l()
        self.dx = self.L / (self.n_x - 1)

        self.x = numpy.linspace(0, self.L, self.n_x)

        self.Y0 = numpy.zeros((4))
        # 2 мерный массив, содержит значения кривизны и её производных в каждой точке
        self.result = None
        # Значение кривизны в каждой точке, одномерный массив
        self.w = numpy.zeros((self.n_x))

    def calculete(self, x, Y):
        """
        Решение дифференциального уравнения с определёнными граничными значениями
        """
        array = numpy.zeros_like(Y)
        array[0] = Y[1]
        array[1] = Y[2]
        array[2] = Y[3]
        # В милиметрах прогиб y
        array[3] = 1000 * self.q / (self.E * self.I)

        return array

    def boundary_condition(self, ya, yb):
        """ Задаёт граничные условия для решения уравнения
        В 0 и на растоянии L координаты и производные равны 0
        """
        return numpy.array([ya[0], yb[0], ya[1], yb[1]])

    def Solve(self):
        '''
        Записывает результат дифференциального уравнения в ячейку объекта self.result
        Записывает значения прогиба в одномерный массив self.w
        '''
        print("Started solving static")
        # Начальное предположение для решения bvp
        guess = numpy.zeros((4, self.x.size))

        res = scipy.integrate.solve_bvp(
            self.calculete, self.boundary_condition, self.x, guess, max_nodes=10000, tol=1e-8)
        self.result = res.sol(self.x)
        self.w = self.result[0]
        print("Ended solving static")

    def Get_render(self):
        """
        Выводит график значений координаты от высоты для каждой точки пластины
        """
        array = numpy.zeros((2 * self.n_x))
        array[:self.n_x] = self.w

        for i in range(self.n_x):  # Пересчёт координаты x
            array[self.n_x + i] = self.dx * i

        matplotlib.pyplot.plot(
            array[self.n_x:], array[:self.n_x], 'b-', linewidth=2)
        matplotlib.pyplot.show()


# Динамическое решение дифференциального уравнения dw2/dt2 = EI/nu dw4/dx4
class Euler_Bernoulli_Dinamic_Finite():
    """
    Решает уравнение колебаний Эйлера - бернулли, методом конечных элементов
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

    def Get_y_x(self, index):
        """
        Выполняет функцию для index значения времени
        Пересчитывает высоту точки через кривизну предыдущех точек
        Возвращает зависимость высоты точки от координаты
        """
        array = numpy.zeros((2 * self.N_x))
        array[:self.N_x] = self.result[index][:self.N_x]  # Y координата

        for i in range(self.N_x):  # Координата x
            array[self.N_x + i] = self.dx * i

        return array

    def Get_render(self, index):
        """
        Выводит визуализацию решения дифференциального уравнения для index значения времени
        """
        array = self.Get_y_x(index)
        y = array[:self.N_x]
        x = array[self.N_x:]
        matplotlib.pyplot.plot(
            x, y, 'b-', linewidth=2)
        matplotlib.pyplot.show()
