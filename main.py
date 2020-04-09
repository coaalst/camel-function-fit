import random
import sys
from decimal import Decimal
import numpy as np
import configparser
import pprint
from matplotlib import pyplot as plt

# Output
#outfile = sys.stdout
outfile = open('output'+'.txt', 'w')

best_ever_sol = None
best_ever_fitment = None


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
def mutiraj(hromo_1, hromo_2, mutation_rate):
    if(random.gauss(0,1) < mutation_rate):
        l = []
        rand = round(random.gauss(-1, 1), 2)
        i = random.randint(0, 1)
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
def function_fit(velicina_populacije, test_vel, broj_pokretanja, best_ever_sol, best_ever_fitment, interval, max_iteracija, mutation_rate, prikazi_graf, save_graf):

    # Podesavanje populacije
    pop_vel = velicina_populacije
    npop_vel = velicina_populacije

    prilagodjenost = []

    best = None
    best_result_fitment = None
    
    for k in range(broj_pokretanja):

        pprint.pprint('Algoritam pokrenut| pokretanje: ' + str(k + 1) + ' od ' + str(broj_pokretanja) + '| populacija: ' + str(velicina_populacije) + '| test velicina: ' + str(test_vel) + "|", outfile)

        # Generisanje populacije pomoću zadatog intervala realnih vrednosti
        pop = [[random.uniform(*interval)
                for i in range(test_vel)] for j in range(pop_vel)]

        i = 0
        t = []
        najbolje_po_gen = []
        prilagodjenost_po_gen = []
        prosek_prilagodjenosti_po_gen = []

        # Ponavljamo dok ne postignemo maksimum iteracija ili dok trošak ne postane 0
        while best_result_fitment != 0 and len(t) != max_iteracija:
            n_pop = pop[:]
            while len(n_pop) < pop_vel + npop_vel:
                h1 = turnir(oceni, pop, 3)
                h2 = turnir(oceni, pop, 3)
                h3, h4 = ukrsti(h1, h2)
                h3, h4 = mutiraj(h3, h4, mutation_rate)

                n_pop.append(h3)
                n_pop.append(h4)

            pop = sorted(n_pop, key=lambda x: oceni(x))[:pop_vel]

            # Dodavanje prilagodjenosti u niz
            for ent in pop:
                prilagodjenost.append(oceni(ent))
                prilagodjenost_po_gen.append(oceni(ent))
            
            prosek = 0
            for l in prilagodjenost_po_gen:
                prosek += l
            prosek_prilagodjenosti_po_gen.append(prosek/len(prilagodjenost_po_gen))
            prilagodjenost_po_gen = []

            # Uzn=mi najboljenje ocenjeni i smesti ga u niz za grafik
            f = oceni(pop[0])
            najbolje_po_gen.append(f)

            # Azuriranje top rezultata
            if best_result_fitment is None or best_result_fitment > f:
                best_result_fitment = f
                best = pop[0]
                # pprint.pprint('Najbolji trenutni: ' + str(best), outfile)
            
            # Append generacije
            t.append(i)
            i += 1
        
        # Grafik
        fig, (ax1, ax2) = plt.subplots(nrows=1, ncols=2)

        # Figure 1
        ax1.plot(t, najbolje_po_gen, color = '#444444')
        ax1.set_xlabel('Generacija')
        ax1.set_ylabel('Najbolji iz generacije')
        
        # Figure 2
        ax2.plot(t, prosek_prilagodjenosti_po_gen)
        ax2.set_xlabel('Generacija')
        ax2.set_ylabel('Prosek prilagodjenosti iz generacije')
        
        if prikazi_graf is True:
            plt.show()

        if save_graf is True: 
            plt.savefig('grafik' + str(best) + '.png', format='png', dpi=500, )


        # Ako smo našli bolji od prethodnog, ažuriramo najbolje rešenje
        if best_ever_fitment is None or best_ever_fitment > best_result_fitment:
            best_ever_fitment = best_result_fitment
            best_ever_sol = best
            pprint.pprint('Najbolji trenutni: ' + str(f), sys.stdout)

        prosecna_prilagodjenost = 0
        for prilag in prilagodjenost:
            prosecna_prilagodjenost += prilag
        
        prosecna_prilagodjenost /= len(prilagodjenost)
        k = 0
        pprint.pprint('Najbolji trenutni fitness: ' + str(oceni(best)), outfile)
        pprint.pprint('Najbolji trenutni: ' + str(best), outfile)
        pprint.pprint('Prosecna prilagodjenost: ' + str(prosecna_prilagodjenost), outfile)


def init():

    # Interval funckije
    interval = []

    # Podesavanje algoritma
    populacije = []
    max_iteracija = 500
    broj_pokretanja = 2
    mutation_rate = 0.2

    test_vel = 0
    pop_vel = 0


    pprint.pprint('Aleksandar Stojanovic RN97/2018', sys.stdout)
    pprint.pprint('Ucitavam konfiguraciju...', sys.stdout)
    config = configparser.ConfigParser()
    
    config_file = input("Unesite ime fajla iz kojeg treba citati parametre, ukoliko zelite default vrednosti upisite (n): ")
    if config_file == 'n':
        config.read_file(open(r'config.properties'))
    else: config.read_file(open(config_file, 'r'))
    
    populacije.append(int(config.get('population', 'size1')))
    populacije.append(int(config.get('population', 'size2')))
    populacije.append(int(config.get('population', 'size3')))

    interval.append(int(config.get('population', 'interval_min')))
    interval.append(int(config.get('population', 'interval_max')))

    mutation_rate = float(config.get('population', 'mutation_rate'))
    test_vel = int(config.get('population', 'test_vel')) 
    max_iteracija = int(config.get('general', 'max_iteracija'))
    broj_pokretanja = int(config.get('general', 'broj_pokretanja'))
    outfile = str(config.get('general', 'output'))
    prikazi_graf = bool(config.get('general', 'prikazi_graf'))
    save_graf = bool(config.get('general', 'save_graf'))

    pprint.pprint('Parser: Ucitana velicina populacije: ' + str(populacije), sys.stdout)
    pprint.pprint('Parser: Ucitan interval: ' + str(interval), sys.stdout)
    pprint.pprint('Parser: Ucitana vrednost koef mutacije: ' + str(mutation_rate), sys.stdout)
    pprint.pprint('Parser: Ucitana velicina testa: ' + str(test_vel), sys.stdout)
    pprint.pprint('Parser: Ucitan broj maks iteracija: ' + str(max_iteracija), sys.stdout)
    pprint.pprint('Parser: Ucitan broj pokretanja: ' + str(broj_pokretanja), sys.stdout)
    pprint.pprint('Parser: Ucitana putanja za output: ' + str(outfile), sys.stdout)

    for velicina_populacije in populacije:
        function_fit(velicina_populacije, test_vel, broj_pokretanja, best_ever_sol, best_ever_fitment, interval, max_iteracija, mutation_rate, prikazi_graf, save_graf)

def main():
    init()


if __name__ == '__main__':
    main()


