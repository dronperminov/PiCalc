#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// метод Валлиса
double vallis_method(double eps) {
	double pi = 2;
	double tmp;
	double a = 2;
	double b = 1;

	do {
		tmp = pi;
		pi *= a * a / (b * (b + 2));
		
		a += 2;
		b += 2;
	} while (fabs(pi - tmp) >= eps);

	return pi;
}

// ряд Нилаканта
double nilakanta_method(double eps) {
	double pi = 3;	
	double tmp;
	double d = 2;
	int sign = 1;

	do {
		tmp = pi;
		pi += sign * 4.0 / (d * (d + 1) * (d + 2));
		
		sign = -sign;
		d += 2;
	} while (fabs(pi - tmp) > eps);

	return pi;
}

// ряд Лейбница
double leibnits_method(double eps) {
	double pi = 0;
	double tmp;
	double n = 0;
	int sign = 1;

	do {
		tmp = pi;
		pi += sign / (2 * n + 1);

		sign = -sign;
		n++;
	} while (fabs(pi - tmp) > eps);

	return pi * 4;
}

// метод Монте-Карло
double monte_carlo_method(double eps) {
	double pi = 0;
	double tmp;

	double total = 32;

	do {
		tmp = pi;
		double within = 0;

		for (long long i = 0; i < total; i++) {
			double x = (double) rand() / RAND_MAX;
			double y = (double) rand() / RAND_MAX;

			if (sqrt(x * x + y * y <= 1))
				within++;
		}

		pi = within / total;

		total *= 2;
	} while (fabs(pi - tmp) > eps);

	return pi * 4;
}

// числовой ряд
double series1_method(double eps) {
	double pi = 0;
	double tmp;
	double n = 0;
	double four = 1;
	int sign = 1;

	do {
		tmp = pi;
		pi += (sign / four) * (2 / (4*n + 1) + 2 / (4*n + 2) + 1 / (4*n + 3));

		sign = -sign;
		four *= 4;
		n++;
	} while (fabs(pi - tmp) > eps);

	return pi;
}

// ещё один числовой ряд
double series2_method(double eps) {
	double pi = 0;
	double tmp;
	double n = 0;
	double three = 1;
	int sign = 1;

	do {
		tmp = pi;
		pi += sign / (three * (2 * n + 1));

		sign = -sign;
		three *= 3;
		n++;
	} while (fabs(pi - tmp) > eps);

	return pi * 2 * sqrt(3);
}

// метод Мэчина
double machin_method(double eps) {
	double pi = 0;
	double tmp;

	double x1 = 1.0 / 5;
	double x2 = 1.0 / 239;
	double y1 = x1;
	double y2 = x2;

	double n = 0;
	int sign = 1;

	do {
		tmp = pi;
		pi += 4 * sign * y1 / (2 * n + 1) - sign * y2 / (2 * n + 1);

		y1 *= x1 * x1;
		y2 *= x2 * x2;
		sign = -sign;
		n++;
	} while (fabs(pi - tmp) > eps);

	return pi * 4;
}

// метод Виета
double viet_method(double eps) {
	double pi = 1;
	double tmp;
	double v = sqrt(2);

	do {
		tmp = pi;
		pi *= v / 2;

		v = sqrt(2 + v);
	} while (fabs(pi - tmp) > eps);

	return 2 / pi;
}

// метод Эйлера через нахождение pi^2/6
double euler_method(double eps) {
	double pi = 0;
	double tmp;
	double n = 1;

	do {
		tmp = pi;
		pi += 1 / (n * n);
		n++;
	} while (fabs(pi - tmp) > eps);

	return sqrt(pi * 6);
}

int main() {
	double eps;
	printf("Enter eps: ");
	scanf("%lf", &eps);

	printf("Method of Nilakanta:   %.16lf\n", nilakanta_method(eps));
	printf("Method of Series 1:    %.16lf\n", series1_method(eps));
	printf("Method of Series 2:    %.16lf\n", series2_method(eps));
	printf("Method of Viet:        %.16lf\n", viet_method(eps));
	printf("Method of Machin:      %.16lf\n", machin_method(eps));
	printf("Method of Leibnits:    %.16lf\n", leibnits_method(eps));
	printf("Method of Vallis:      %.16lf\n", vallis_method(eps));
	printf("Method of Euler:       %.16lf\n", euler_method(eps));
	printf("Method of Monte Carlo: %.16lf\n", monte_carlo_method(eps));
}