import random
import sys


def predvidi(hromozom, x):
      return hromozom[0] ** x + hromozom[1]

# Ocena gena
def trosak(hromozom):
  greska = 0
  for x, y in enumerate(podaci, 1):
    predvidjanje = predvidi(hromozom, x)
    greska += (predvidjanje - y[1]) ** 2
  return greska / len(podaci)

# Funckija za mutiranje - Tackasta normalna mutacija za kontinualni
# genetrski algoritam (dodaje se slucajna vrednost na gen iz normalne raspodele).
def mutiraj(hromozom, rate, opseg):
    if random.random() < rate:
        for i in range(len(hromozom)):
#            hromozom[i] = random.uniform(*opseg)
            hromozom[i] += random.gauss(0, 1)
    return hromozom

#p = p + sigmaT(0,1)

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

def function_fit():

    #Podesavanje populacije
    pop_vel = 20
    npop_vel = 20
    max_iteracija = 500
    mutacije = [0.05, 0.1, 0.2]

    #interval funckije
    interval = [-5, 5]
    
    outfile = sys.stdout
    s_trosak = 0
    s_iteracija = 0
    best_ever_sol = None
    best_ever_f = None
    test_vel = 2
    # 5 pokretanja genetskog algoritma
    for k in range(max_iteracija):
        print('Algoritam pokrenut, iteracija: ' + k, ', mutacija: ' + mut_rat, file=outfile)
        best = None
        best_f = None
        t = 0

        # generisanje populacije pomoću zadatog intervala realnih vrednosti
        pop = [[random.uniform(*interval) for i in range(test_vel)] for j in range(pop_vel)]
        # ponavljamo dok ne postignemo maksimum iteracija ili dok trošak ne postane 0
        while best_f != 0 and t < max_iteracija:
            n_pop = pop[:]
            while len(n_pop) < pop_vel + npop_vel:
                h1 = turnir(trosak, pop, 3)
                h2 = turnir(trosak, pop, 3)
                h3, h4 = ukrsti(h1, h2)
                mutiraj(h3, mut_rat, opseg)
                mutiraj(h4, mut_rat, opseg)
                n_pop.append(h3)
                n_pop.append(h4)
            pop = sorted(n_pop, key=lambda x : trosak(x))[:pop_vel]
            f = trosak(pop[0])
            if best_f is None or best_f > f:
                best_f = f
                best = pop[0]
#                    print(t, best_f, file=outfile)
            t += 1
        s_trosak += best_f
        s_iteracija += t
        # ako smo našli bolji od prethodnog, ažuriramo najbolje rešenje
        if best_ever_f is None or best_ever_f > best_f:
            best_ever_f = best_f
            best_ever_sol = best
        print(t, best, best_f, file=outfile)
#                print(t, best, best_f)
    # na kraju svih izvršavanja izračunavamo srednji trošak i srednji broj iteracija
    s_trosak /= 5
    s_iteracija /= 5
    print('Srednji trosak: %g' % s_trosak, file=outfile)
    print('Srednji broj iteracija: %.2f' % s_iteracija, file=outfile)
    print('Najbolje resenje: %s' % best_ever_sol, file=outfile)
    print('Najbolji trosak: %g' % best_ever_f, file=outfile)
    print('Procena sledećeg unosa: %d' % int(predvidi(best_ever_sol, len(podaci) + 1)))


for mut_rat in mutacije:
    function_fit()