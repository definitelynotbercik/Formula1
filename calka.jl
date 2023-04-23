using Polynomials


####### Przykładowe trasy ######
using Distributions
using Plots

ograniczenie_d = 10 #promien krótszego łuku
ograniczenie_g = 20 #promien dluzszego luku
punktp_sx = 15 #z tego punktu start (x = 15)
punktp_sy = 0 #z tego punktu start (y =0)

#tworzenie trasy
f(x) = ograniczenie_g*sin(x)
g(x) = ograniczenie_g*cos(x)
h(x) = ograniczenie_d*sin(x)
j(x) = ograniczenie_d*cos(x)
tab1 = [g,j]
tab2 = [f,h]

tab3 = Array{Function}(undef, 0)
tab4 =Array{Function}(undef, 0)

x3 =[]
y  =[]
r  =[]
al = []
dl =[]
k = rand(Uniform(-300,0), 1000) #1000 przypadkowych wartosci z przedzialu;
#wazne zeby drugie bylo 0, bo wychodza dziwne zakrety
c = rand(Uniform(-300,300), 1000)
for t1 in k
    for t2 in c
        x2 = t1
        y3 =x2 + t2
        r1 = sqrt((x2-punktp_sx)^2 + (y3)^2) #promien danego okregu
        odl = sqrt((x2)^2+ (y3)^2) #odl miedzy punktem (0,0), a srodkiem danego okregu
        y0 = y3 + sqrt((r1)^2 + (x2)^2) # punkt przeciecia z osia OY (0, y0)
        if( abs(r1-ograniczenie_d)>odl && y0<ograniczenie_g) #ograniczenia
            al1 = acos(sqrt(7.5^2+(y0/2)^2)/r1) #kąt tworzony przez punkty: (0,y0), sr danego okregu i (15,0) 
            dl1 = al1*r1  #dlugosc luku okregu miedzy (0,y0) i (15,0) 
            push!(y, y3)
            push!(x3, x2)
            push!(r, r1)
            push!(al, al1)
            push!(dl, dl1)
        end
    end
end

for s in 1:length(r)
    z(x) = r[s]*sin(x) + y[s]
    w(x) = r[s]*cos(x) + x3[s]
    push!(tab3, w)
    push!(tab4, z)
end

print(argmin(dl), "\n") #numer argumentu o najkrotszym luku
print(dl[argmin(dl)], "\n") #dlugosc najkrotszego luku
plot(tab3,tab4,-pi,pi, aspectratio=1, linewidth=1, xlim = (0, ograniczenie_g),
ylim = (0, ograniczenie_g), legend = false)
plot!(tab1,tab2,0,pi/2, aspectratio=1, linewidth=5, thickness_scaling = 1, c=:black )


function wielomian(r)
    # Wyznaczenie punktów końcowych łuku
    x1, y1 = r, 0
    x2, y2 = 0, r

    # Wyznaczenie równania prostej przechodzącej przez punkty (x1, y1) i (x2, y2)
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1

    # Wyznaczenie punktów przecięcia okręgu i prostej
    a = m^2 + 1
    b = 2 * m * b
    c = b^2 / 4 - r^2
    delta = sqrt(b^2 - 4 * a * c)
    x_1 = (-b + delta) / (2 * a)
    #x_2 = (-b - delta) / (2 * a)
    y_1 = m * x_1 + b
    #y_2 = m * x_2 + b

    # Wyznaczenie długości łuku między punktem (x1, y1) a punktem przecięcia (x_1, y_1)
    d = sqrt((x_1 - x1)^2 + (y_1 - y1)^2)

    # Wyznaczenie liczby punktów na łuku
    n = Int(round((d / r) * pi / 2))

    # Wyznaczenie punktów na łuku
    pary_punktow_luku = [(r * cos(t), r * sin(t)) for t in range(0, stop=pi/2, length=n)]

    # Wyznaczenie stopnia wielomianu
    stopien = div(n, 2)

    # Wyznaczenie współczynników wielomianu
    xs = [point[1] for point in pary_punktow_luku]
    ys = [point[2] for point in pary_punktow_luku]
    A = ones(n, stopien+1)
    for i=1:n
        for j=2:(stopien+1)
            A[i,j] = A[i,j-1] * xs[i]
        end
    end
    punkty = (A' * A) \ (A' * ys)
    # Zwrócenie wielomianu
    return punkty
end

function calculate_distance(f::Function)
    f_prime = ForwardDiff.derivative(f)    #Pochodna wielomianu określającego 
    return sqrt(1 + f_prime(x)^2)
end

wiel=(Polynomial(wielomian(r[argmin(dl)]), :x))
print(wiel, "\n")
r_nowe=[]
for i in r
   push!(r_nowe, wielomian(i))
end
print(derivative(wiel), "\n")
pochodna=((((derivative(wiel))^2)+1)) #brakuje pierwiastka ale nie wiem jak to zrobić
calka = integrate(pochodna)
print(calka, "\n")
print(calka(r[argmin(dl)])-calka(0))