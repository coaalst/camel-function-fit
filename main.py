import random
import sys
from decimal import Decimal
import numpy as np

# Ocena gena se vrsi tako sto za svaki hromozom(najbolji) gledamo koliko je blizu dosao
# ocekivanom rezultatu.
# 5.3 Funkcija prilagodjenosti
# Уколико је проблем оптимизације такав да је потребно оптимизовати реалну функцију више променљивих,
# онда се она сама може искористити као функција прилагођености. У том случају заобилази се и корак
# моделовања функције прилагођености, док се њени параметри директно преписују у хромозом. На тај
# начин употреба генетског алгоритма за оптимизацију вишедимензионалне реалне функције
# постаје значајно олакшана.
def oceni(hromozom):

    # Ocekivani minimum f(0, 0) = 0
    x = float(hromozom[0])
    y = float(hromozom[1])
    return 2 * pow(x, 2) - 1.05 * pow(x, 4) + (pow(x, 6) / 6) + x * y + y * y


# Funckija za mutiranje - Tackasta normalna mutacija za kontinualni
# genetrski algoritam (dodaje se slucajna vrednost na gen iz normalne raspodele).
def mutiraj(hromo_1, hromo_2, rate, opseg ):
    if random.random() < rate:
        l = []
        rand = round(random.gauss(0,1), 2)
        i = random.randint(0,1)
        hromo_1[i] = float(hromo_1[i]) + rand
        hromo_2[i] = float(hromo_2[i]) + rand
        l.append(hromo_1)
        l.append(hromo_2)
        return l
    return hromo_1, hromo_2


# turnirska selekcija - argumenti su funkcija troška, rešenje, populacija i veličina turnira
def turnir(fja, pop, vel):
    z = []
    while len(z) < vel:
        z.append(random.choice(pop))
    najbolji = None
    najbolji_f = None
    for e in z:
        ff = fja(e)
        if najbolji is None or ff < najbolji_f:
            najbolji_f = ff
            najbolji = e
    return najbolji


# Metoda ukrstanja - jednostavno ukrstanje REDKLIF
# Редклиф је у свом раду предложио једноставан облик формуле за израчунавање вредности једног
# параметра потомка из вредности истог параметра његова два родитеља [19]:
# 𝑝pot=𝛽𝑝1+(1−𝛽)𝑝2
# У овој формули 𝑝pot представља параметар потомка, 𝑝1 i 𝑝2 су параметри првог и другог родитеља
# респективно, а 𝛽∈[0,1] је случајни број којим се одређује удео вредности параметра првог и
# другог родитеља у вредности параметра добијеног потомка. Када је 𝛽=0, вредност параметра се
# преписује од другог родитеља; када је 𝛽=1, вредност се преписује из првог родитеља, а када је
# 𝛽=0.5, вредност је аритметичка средина вредности параметра оба родитеља.
def ukrsti(hromo_1, hromo_2):
    beta_fact = []
    beta_fact.append(random.gauss(0, 1))
    beta_fact.append(random.gauss(0, 1))
    l = [[], []]
    for i in range(len(l)):
        x = beta_fact[i] * hromo_1[0] + ((1 - beta_fact[i]) * hromo_2[0])
        y = beta_fact[i] * hromo_1[1] + ((1 - beta_fact[i]) * hromo_2[1])
        l[i] = [x, y] 
    return l

# Glavna funckija koja kao argument uzima svaki rejt mutacije iz niza
def function_fit(mut_rat):

    # Podesavanje populacije
    pop_vel = 20
    npop_vel = 20
    max_iteracija = 500

    # Interval funckije
    interval = [-5, 5]

    # Output
    outfile = sys.stdout

    # Podesavanje algoritma
    s_oceni = 0
    s_iteracija = 0
    best_ever_sol = None
    best_ever_f = None
    test_vel = 2

    # Pokrecemo 2 puta po mutaciji
    for k in range(2):
        print('Algoritam pokrenut, iteracija: ' + str(k),
              ', mutacija: ' + str(mut_rat), file=outfile)
        best = None

        # Generisanje populacije pomoću zadatog intervala realnih vrednosti
        pop = [[round(random.uniform(*interval),2) for i in range(test_vel)] for j in range(pop_vel)]
        # print('x =' + str(pop.pop()[0]) + ', y =' + str(pop.pop()[1]))

        best_result_fitment = None
        t = 0
        # Ponavljamo dok ne postignemo maksimum iteracija ili dok trošak ne postane 0
        while best_result_fitment != 0 and t < max_iteracija:
            n_pop = pop[:]
            while len(n_pop) < pop_vel + npop_vel:
                h1 = turnir(oceni, pop, 3)
                h2 = turnir(oceni, pop, 3)
                h3, h4 = ukrsti(h1, h2)
                h3, h4 = mutiraj(h3, h4, mut_rat, interval)

                n_pop.append(h3)
                n_pop.append(h4)

            pop = sorted(n_pop, key=lambda x: oceni(x))[:pop_vel]
            f = oceni(pop[0])
            if best_result_fitment is None or best_result_fitment > f:
                best_result_fitment = f
                best = pop[0]
                print('Najbolji trenutni: ' + str(best))
            t += 1

        # Azuriraj global statistiku
        s_oceni += best_result_fitment
        s_iteracija += t

        # Ako smo našli bolji od prethodnog, ažuriramo najbolje rešenje
        if best_ever_f is None or best_ever_f > best_result_fitment:
            best_ever_f = best_result_fitment
            best_ever_sol = best
        print('Mutacija:', t, best, best_result_fitment, file=outfile)

    # Na kraju svih izvršavanja izračunavamo srednji broj iteracija
    s_iteracija /= 2
    k = 0
    print('Srednji broj iteracija: %.2f' % s_iteracija, file=outfile)
    print('Najbolje resenje: %s' % best_ever_sol, file=outfile)


mutacije = [0.05, 0.1, 0.2]

print('Aleksandar Stojanovic RN97/2018')
for mut_rat in mutacije:
    function_fit(mut_rat)
