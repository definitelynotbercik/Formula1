#wbudowany pakiet
#metoda trapezów
using QuadGK
function integrand(s::Float64, r1::Float64, r2::Float64, s1::Float64, m::Float64, g::Float64)
    u=r1+(r2-r1)/s1*s
    sqrt(u/m*g)  
end

function integrate(r1::Float64, r2::Float64, s1::Float64, m::Float64, g::Float64, a::Float64, b::Float64)
    f(s) = integrand(s, r1, r2, s1, m, g)
    Q, err=quadgk(f, a, b)
    Q, err
end
i = (integrate(2.0, 7.0, 4.0, 9.0, 2.0, 2.0, 10.0))

#funkcja ręcznie pisana
#metoda prostokątów
function integrand(s::Float64, r1::Float64, r2::Float64, s1::Float64, m::Float64, g::Float64)
    v = r1 + (r2 - r1) / s1 * s
    return sqrt(v / m * g)
end

function oblicz_calka(a::Float64, b::Float64)
    n = 1000
    dx = (b - a) / n
    x = a + dx / 2 : dx : b - dx / 2
    Q = sum(integrand.(x, 2.0, 7.0, 4.0, 9.0, 2.0)) * dx
    return Q
end

i = (oblicz_calka(2.0, 10.0))




# ##GENEROWANIE ZAKRĘTU
# using Plots

# function zakret()
#     r_min=randrange(0, 10)
#     r_max=randrange(0, 10)

#     theta = range(0, stop=2π, length=n_points)
#     r = rand(n_points) * (r_max - r_min) + r_min
#     return [(r[i] * cos(theta[i]), r[i] * sin(theta[i])) for i in 1:n_points]
# end
# curve = zakret()

# function fit_circle_to_curve(curve::Vector{Tuple{Float64, Float64}})::Tuple{Float64, Float64, Float64}
#     x = [p[1] for p in curve]
#     y = [p[2] for p in curve]
#     x_mean = mean(x)
#     y_mean = mean(y)
#     u = x .- x_mean
#     v = y .- y_mean
#     Suu = sum(u.^2)
#     Svv = sum(v.^2)
#     Suv = sum(u .* v)
#     Suvv = sum(u .* v.^2)
#     Svuu = sum(v .* u.^2)
#     a = (Svv * (Suu + Svv) - Suv^2) / (2 * Svv * u' * u - 2 * Suv * u' * v)
#     b = (Suu * (Suu + Svv) - Suv^2) / (2 * Suu * v' * v - 2 * Suv * u' * v)
#     r = sqrt(a^2 + b^2 + (Suu + Svv) / length(x))
#     return (x_mean, y_mean, r)
# end

# function plot_curve_and_circle(curve::Vector{Tuple{Float64, Float64}}, circle::Tuple{Float64, Float64, Float64})
#     plot(curve, seriestype=:scatter)
#     x_mean, y_mean, r = circle
#     x = x_mean .+ r * cos.(range(0, stop=2π, length=100))
#     y = y_mean .+ r * sin.(range(0, stop=2π, length=100))
#     plot!(x, y)
# end

# curve = zakret()
# circle = fit_circle_to_curve(curve)

##GENETROWANIE ZAKRĘTU O DANYM PROMIENIU
using LinearAlgebra
using Plots
gr()

r = 5.0  # promień łuku
r2=3.0
arc_length = π/2 # długość łuku
num_points = 100 # liczba punktów

theta = range(0.0, arc_length, num_points) # równomierne wartości kąta
x = r * cos.(theta)
y = r * sin.(theta)
x2 = r * cos.(theta)/2
y2 = r * sin.(theta)/2


p1=plot(x, y, seriestype =:line, aspect_ratio = 1)
p2=plot(x2, y2, seriestype =:line, aspect_ratio = 1)
plot([x x2], [y y2])


# ## Tor jako 1/4 elipsy
# using Plots

# # Długości półosi elipsy
# a = 3
# b = 0.5
# d = 2
# e = 1
# # Centrum elipsy
# c = -0.7 + 0.27im

# # Funkcja zwracająca punkt (x, y) na elipsie dla kąta theta
# get_point(theta) = (real(c) + a*cos(theta), imag(c) + b*sin(theta))
# get_point(theta2) = (real(c) + d*cos(theta2), imag(c) + e*sin(theta2))
# # Generowanie punktów na elipsie dla kątów od 0 do pi/2
# thetas = range(0, stop=pi/2, length=100)
# points = [get_point(theta) for theta in thetas]
# points2 = [get_point(theta2) for theta2 in thetas]
# # Rysowanie wykresu
# p1=plot([x for (x,y) in points], [y for (x,y) in points], legend=false, aspect_ratio=:equal, linewidth=2)
# p2=plot([x for (x,y) in points2], [y for (x,y) in points2], legend=false, aspect_ratio=:equal, linewidth=2)
# plot(p1,p2, legend=false)

using LinearAlgebra
using Plots
using Statistics
gr()

r = 5.0  # promień łuku
r2 = 3.0
arc_length = π/2 # długość łuku
num_points = 100 # liczba punktów

theta = range(0.0, arc_length, num_points) # równomierne wartości kąta
x = r * cos.(theta)
y = r * sin.(theta)
x2 = r2 * cos.(theta)
y2 = r2 * sin.(theta)

# generowanie losowych wartości dla każdego punktu
randvalue = rand(num_points)
xrand=[]
yrand=[]
# skalowanie losowych wartości przez długość łuku i dodawanie do punktów na łuku
for i in 1:2
    x_random = x.+randvalue*(x2[i]-x[i])
    y_random = y.+randvalue*(y2[i]-y[i])
    push!(xrand, x_random)
    push!(yrand, y_random)
end

# # rysowanie łuków i losowej trasy
# plot(x, y, seriestype = :line, aspect_ratio = 1, label = "łuk 1")
# plot!(x2, y2, seriestype = :line, aspect_ratio = 1, label = "łuk 2")
# plot!(xrand, yrand, seriestype = :line, label = "losowa trasa")




using Polynomials

# generowanie losowych wartości dla każdego punktu
randvalue = rand(num_points)
xrand = []
yrand = []

# skalowanie losowych wartości przez długość łuku i dodawanie do punktów na łuku
for i in 1:2
    x_random = x .+ randvalue .* (x2[i] - x[i])
    y_random = y .+ randvalue .* (y2[i] - y[i])
    push!(xrand, x_random)
    push!(yrand, y_random)
end

# rysowanie łuków i losowej trasy
plot(x, y, seriestype = :line, aspect_ratio = 1, label = "łuk 1")
plot!(x2, y2, seriestype = :line, aspect_ratio = 1, label = "łuk 2")
plot!(xrand, yrand, seriestype = :line, label = "losowa trasa")

# wyznaczanie trendline
degree = 1  # stopień wielomianu
coeffs = fit(xrand[1], yrand[1], degree)  # współczynniki wielomianu dla pierwszej trajektorii
x_trend = range(minimum(xrand[1]), maximum(xrand[1]), length=100)
y_trend = evalpoly(coeffs, x_trend)

plot(x_trend, y_trend, linewidth=2, linestyle=:dash, label="trendline")



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
        dluod = sqrt(15^2+y0^2) #odl miedzy punktem (15,0) a (0,y0)
        if( abs(r1-ograniczenie_d)>odl && y0<ograniczenie_g) #ograniczenia
            al1 = acos(1-((dluod^2)/(2*r1^2))) #kąt tworzony przez punkty: (0,y0), sr danego okregu i (15,0) 
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
print(dl[argmin(dl)]) #dlugosc najkrotszego luku
plot(tab3,tab4,-pi,pi, aspectratio=1, linewidth=1, xlim = (0, ograniczenie_g),
ylim = (0, ograniczenie_g), legend = false)
plot!(tab1,tab2,0,pi/2, aspectratio=1, linewidth=5, thickness_scaling = 1, c=:black )


#### WSTEPNA ANIMACJA DLA NAJPROSTSZEJ TRASY ####
using Plots

@userplot CirclePlot
@recipe function f(cp::CirclePlot)
    x, y, i = cp.args
    n = length(x)
    inds = circshift(1:n, 1 - i)
    linewidth --> range(0, 10, length = n)
    seriesalpha --> range(0, 1, length = n)
    aspect_ratio --> 1
    label --> false
    x[inds], y[inds]
end

n = 150
t = range(0, 2π, length = n)
x = 15*sin.(-t)
y = 15*cos.(-t)

anim = @animate for i ∈ 1:n
    circleplot!(x, y, i,xlim=(0,1), ylim=(0, 1))
end
@gif for i ∈ 1:n
    circleplot(x, y, i, line_z = 1:n, cbar = false, framestyle = :zerolines, xlim=(0,20), ylim=(0, 20))
    plot!([10*cos.(-t), 20*cos.(-t)], [10*sin.(-t), 20*sin.(-t)], )
end every 5
