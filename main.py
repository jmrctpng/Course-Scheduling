from typing import List
import random
import time
import numpy as np
import copy
import threading

WH_CACHE = []
PL_CACHE = []
BL_CACHE = []
PWL_CACHE = []
MEMOIZATION = True


class Data:
    professor = []
    block = []
    lab_room = []
    reg_room = []
    lab_course = []
    lec_course = []

    @staticmethod
    def add_lab_course(course):
        Data.lab_course.append(course)

    @staticmethod
    def add_lec_course(course):
        Data.lec_course.append(course)

    @staticmethod
    def add_room(room):
        if room.type_of_room == "Regular":
            Data.reg_room.append(room)
        else:
            Data.lab_room.append(room)

    @staticmethod
    def add_prof(prof):
        Data.professor.append(prof)

    @staticmethod
    def add_block(block):
        Data.block.append(block)

    @staticmethod
    def find_prof(professor):
        for prof in Data.professor:
            if prof.get_prof_id() == professor.get_prof_id():
                return prof
        return None

    @staticmethod
    def find_room(room_id, type_of_room):
        rooms = Data.reg_room if type_of_room == "Regular" else Data.lab_room
        for room in rooms:
            if room.get_code() == room_id:
                return room
        return None

    @staticmethod
    def find_block(block):

        for b in Data.block:
            if b.get_block_code() == block.get_block_code():
                return block

        return None

    @staticmethod
    def find_course(course_id):

        course_list = []
        course_list.extend(Data.lab_course)
        course_list.extend(Data.lec_course)

        for course in course_list:
            if course.get_code() == course_id:
                return course

        return None

    @staticmethod
    def find_lec_course(course_id):

        for course in Data.lec_course:
            if course.get_code() == course_id:
                return course

        return None

    @staticmethod
    def find_lab_course(course_id):

        for course in Data.lab_course:
            if course.get_code() == course_id:
                return course

        return None


class Schedule:
    def __init__(self):
        self.schedule = {
            'M': [],
            'T': [],
            'W': [],
            'Th': [],
            'F': [],
            'S': []
        }

    def set_schedule(self, day, time):
        sched = self.schedule[day]
        sched.append(time)

    def clear_and_add_new_sched(self, day, sched):
        print("clear : ", sched)
        print("Before update:", self.schedule[day])
        self.schedule[day] = copy.copy(sched)
        print("After update:", self.schedule[day])

    def get_schedule(self, day):
        return self.schedule[day]

    def remove_sched(self, day, value):
        print("x value : ", value)
        print("x sched : ", self.schedule[day])
        self.schedule[day].remove(value)
        print("x sched11 : ", self.schedule[day])

    def get_all_schedules(self):
        return self.schedule

    def get_m(self):
        return self.schedule['M']

    def get_t(self):
        return self.schedule['T']

    def get_w(self):
        return self.schedule['W']

    def get_th(self):
        return self.schedule['Th']

    def get_f(self):
        return self.schedule['F']

    def get_s(self):
        return self.schedule['S']


class GeneticAlgoPopulation:
    chromosome = []


class GeneticAlgorithm:

    def __init__(self):
        self.chromosome = Chromosome()
        self.new_chromosome = Chromosome()
        self.chromosome_size = 0
        self.curr_chromosome = None
        self.class_slot_population = None

    def encode_chromosome(self, population_class):
        """
        Randomly select one item from each population list and create a new chromosome
        list containing those four items. The for loop iterates over the length of
        self.course_population and creates one chromosome for each course.
        """

        data = Data()

        print("CS size : ", len(population_class.get_class_slots()))
        print("course size : ", len(population_class.get_courses()))
        print("block size : ", len(population_class.get_blocks()))
        print("prof size : ", len(population_class.get_professors()))

        for i in range(len(population_class.get_courses())):
            class_slot = population_class.get_class_slots().pop(
                random.randrange(len(population_class.get_class_slots())))
            course = population_class.get_courses().pop(random.randrange(len(population_class.get_courses())))
            block = population_class.get_blocks().pop(random.randrange(len(population_class.get_blocks())))
            professor = population_class.get_professors().pop(random.randrange(len(population_class.get_professors())))

            gene = Gene(class_slot, course, block, professor)

            # set schedule for professor
            found_prof = data.find_prof(professor)
            if found_prof is not None:
                found_prof.set_schedule(class_slot[1], class_slot[2])
            else:
                print("prof not found")

            # set schedule for classroom
            found_room = data.find_room(class_slot[0].get_code(), class_slot[0].get_type_of_room())
            if found_room is not None:
                found_room.set_schedule(class_slot[1], class_slot[2])

            # set schedule for block
            found_block = data.find_block(block)
            if found_block is not None:
                found_block.set_schedule(class_slot[1], class_slot[2])

            self.chromosome.add_gene(gene)

        self.class_slot_population = population_class.get_class_slots()
        GeneticAlgoPopulation.chromosome.append(self.chromosome)

    def calculate_fit(self):
        # get the current chromosome

        print("Calculate fit size : ", len(GeneticAlgoPopulation.chromosome))
        if len(GeneticAlgoPopulation.chromosome) == 1:
            self.curr_chromosome = GeneticAlgoPopulation.chromosome[0]
        else:
            self.curr_chromosome = GeneticAlgoPopulation.chromosome[-1]

        sum = 0
        self.chromosome_size = len(self.curr_chromosome.get_genes())

        print("X Chromosome Size : ", self.chromosome_size)

        print("Current Chromosome Genes : ", self.curr_chromosome.get_genes())

        # calculate gene fitness score
        for gene in self.curr_chromosome.get_genes():
            gene_fitness = GeneFitness(gene)
            fitness_score = self.fitness_value(gene_fitness)

            if gene_fitness.first_year() is None:
                first_year = 1
            else:
                first_year = gene_fitness.first_year()

            gene.add_fitness_score(fitness_score)
            sum = sum + fitness_score

            print("Gene ID : ", gene)
            print("Class slot: ", gene.class_slot[0].get_code(), " ", gene.class_slot[0].get_capacity(), " ",
                  gene.class_slot[0].get_type_of_room(), " ", gene.class_slot[1], " ", gene.class_slot[2])
            print("Course code: ", gene.course.get_code())
            print("Course hr : ", gene.course.get_hour())
            print("Professor: ", gene.professor.get_prof_id())
            print("Block: ", gene.block.get_block_code())
            print("Score :", gene.fitness_score)
            print("CR :", gene.get_classroom_capacity_fitness())
            print("PW :", gene.get_professor_work_load_fitness())
            print("RA :", gene.get_room_availability_fitness())
            print("BA :", gene.get_block_availability_fitness())
            print("RS :", gene.get_room_suitability_fitness())
            print("CSS :", gene.get_course_slot_suitability_fitness())
            print("FF :", first_year)
            print("PL :", gene.get_prof_lunch_break_fitness())
            print("BL :", gene.get_block_lunch_break_fitness())
            print("WH :", gene.get_working_hours_fitness())
            print("PHC :", gene.get_prof_handled_course_fitness())
            print("BEC :", gene.get_block_enrolled_course_fitness())
            print("=======================================")
            print("prof id : ", gene.professor)
            print("prof sched : ", gene.professor.get_schedule().get_schedule(gene.class_slot[1]))
            print("block sched : ", gene.block.get_schedule().get_schedule(gene.class_slot[1]))
            print("Room sched : ", gene.class_slot[0].get_schedule().get_schedule(gene.class_slot[1]))
            print("-----------------------------")

        print("summm : ", sum, " get_len_gene : ", self.chromosome_size)
        print("new chrooom : ", len(self.new_chromosome.get_genes()))
        self.curr_chromosome.set_fitness_value(sum / self.chromosome_size)

        print("chromosome total fit : ", self.curr_chromosome.get_fitness_value())

        return self.curr_chromosome.get_fitness_value()

    def roulette_wheel(self):

        print("Current Chromosome Size : ", len(self.curr_chromosome.get_genes()))

        index = 0
        selected_gene = None

        if self.chromosome_size > 0:

            random_number = np.random.random()

            while random_number > 0:
                if index < self.chromosome_size:
                    gene = self.curr_chromosome.get_specific_gene(index)
                    gene_score = gene.get_fitness_score()
                    random_number = random_number - gene_score
                    index = index + 1
                else:
                    break

            index = index - 1
            selected_gene = self.curr_chromosome.get_genes().pop(index)
            # self.curr_chromosome.pop_gene(index)
            self.chromosome_size = self.chromosome_size - 1
            print("update size : ", self.chromosome_size)

        return selected_gene

    def hill_climbing(self, parent1):

        attributes_to_improve = GeneticAlgorithm.check_attribute(parent1)

        if len(attributes_to_improve) == 0:
            return None
        else:
            if self.chromosome_size > 0:

                best = self.parent_two(self.curr_chromosome.get_genes(), attributes_to_improve, parent1)
                return best
            else:
                return 0

    @staticmethod
    def check_attribute(parent1):

        search_attribute = []

        if parent1.get_professor_work_load_fitness() == 0 or parent1.get_room_availability_fitness() == 0 or parent1.get_block_availability_fitness() == 0 or parent1.get_course_slot_suitability_fitness() == 0 or parent1.get_prof_lunch_break_fitness() == 0 or parent1.get_block_lunch_break_fitness() == 0:
            search_attribute.append("cs")  # change course slot
        if (
                parent1.get_prof_handled_course_fitness() == 1 and parent1.get_block_enrolled_course_fitness() == 1 and parent1.get_room_suitability_fitness() == 0) or (
                parent1.get_prof_handled_course_fitness() == 0 and parent1.get_block_enrolled_course_fitness() == 0 and parent1.get_room_suitability_fitness() == 1) or (
                parent1.get_classroom_capacity_fitness() == 0):
            search_attribute.append("r")  # change room
        if (
                parent1.get_prof_handled_course_fitness() == 0 and parent1.get_block_enrolled_course_fitness() == 0 and parent1.get_room_suitability_fitness() == 0) or (
                parent1.get_prof_handled_course_fitness() == 1 and parent1.get_block_enrolled_course_fitness() == 0 and parent1.get_room_suitability_fitness() == 0):
            search_attribute.append("c")  # change course
        if parent1.get_prof_handled_course_fitness() == 0 and parent1.get_block_enrolled_course_fitness() == 1 and parent1.get_room_suitability_fitness() == 1:
            search_attribute.append("p")  # change prof
        if (
                parent1.get_prof_handled_course_fitness() == 1 and parent1.get_block_enrolled_course_fitness() == 0 and parent1.get_room_suitability_fitness() == 1) or (
                parent1.get_prof_handled_course_fitness() == 0 and parent1.get_block_enrolled_course_fitness() == 1 and parent1.get_room_suitability_fitness() == 0):
            search_attribute.append("b")  # change block
        if parent1.get_first_year_fitness() == 0 or parent1.get_working_hours_fitness() == 0:
            search_attribute.append("d")  # change day

        return search_attribute

    def parent_two(self, items, attribute, parent1):

        perfect_gene_found = None
        perfect_gene_score = 0
        prev_time_slot = parent1.get_class_slot_one(2)
        prev_day = parent1.get_class_slot_one(1)

        orig_block_sched = copy.copy(parent1.get_block().get_schedule().get_schedule(prev_day))
        orig_prof_sched = copy.copy(parent1.get_professor().get_schedule().get_schedule(prev_day))

        print("P2: P1 ID : ", parent1)
        print("P2: P1 TS : ", prev_time_slot)
        print("P2: P1 BS : ", parent1.get_block().get_schedule().get_schedule(prev_day))
        print("P2: P1 PS : ", parent1.get_professor().get_schedule().get_schedule(prev_day))
        print("P2: P1 CS : ", parent1.get_class_slot_one(0).get_schedule().get_schedule(prev_day))

        for i in range(12):
            random_num = np.random.randint(self.chromosome_size)

            item = items[random_num]

            different_day = False
            remove_prof_sched = False
            remove_block_sched = False

            parent1_copy = copy.copy(parent1)
            item_prev_day = item.get_class_slot_one(1)

            item_orig_prof_sched1 = copy.copy(item.get_professor().get_schedule().get_schedule(item_prev_day))
            item_orig_block_sched1 = copy.copy(item.get_block().get_schedule().get_schedule(item_prev_day))
            item_orig_prof_sched2 = None
            item_orig_block_sched2 = None
            orig_block_sched2 = None
            orig_prof_sched2 = None

            if prev_day != item_prev_day:
                item_orig_prof_sched2 = copy.copy(item.get_professor().get_schedule().get_schedule(prev_day))
                item_orig_block_sched2 = copy.copy(item.get_block().get_schedule().get_schedule(prev_day))
                orig_block_sched2 = copy.copy(parent1.get_block().get_schedule().get_schedule(item_prev_day))
                orig_prof_sched2 = copy.copy(parent1.get_professor().get_schedule().get_schedule(item_prev_day))
                different_day = True

            if "cs" in attribute:
                parent1_copy.set_class_slot(item.get_class_slot())

                new_time_slot = parent1_copy.get_class_slot_one(2)
                new_day = parent1_copy.get_class_slot_one(1)

                if prev_day == new_day and prev_time_slot == new_time_slot:
                    pass
                else:
                    parent1_copy.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                    parent1_copy.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)

                    remove_prof_sched = True
                    remove_block_sched = True

                    parent1_copy.get_block().set_schedule(new_day, new_time_slot)
                    parent1_copy.get_professor().set_schedule(new_day, new_time_slot)

            if "r" in attribute:
                if item.get_class_slot_one(0).get_code() != parent1.get_class_slot_one(0).get_code():
                    parent1_copy.set_class_slot(item.get_class_slot())

                    new_time_slot = parent1_copy.get_class_slot_one(2)
                    new_day = parent1_copy.get_class_slot_one(1)

                    if prev_day == new_day and prev_time_slot == new_time_slot:
                        pass
                    else:

                        if not remove_block_sched:
                            parent1.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                            parent1_copy.get_block().set_schedule(new_day, new_time_slot)
                            remove_block_sched = True
                        if not remove_prof_sched:
                            parent1.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                            parent1_copy.get_professor().set_schedule(new_day, new_time_slot)
                            remove_prof_sched = True

            if "d" in attribute:
                if item.get_class_slot_one(1) != parent1.get_class_slot_one(1):
                    parent1_copy.set_class_slot(item.get_class_slot())

                    new_time_slot = parent1_copy.get_class_slot_one(2)
                    new_day = parent1_copy.get_class_slot_one(1)

                    # remove prev sched and append new sched
                    if not remove_block_sched:
                        parent1.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                        parent1_copy.get_block().set_schedule(new_day, new_time_slot)
                        remove_block_sched = True
                    if not remove_prof_sched:
                        parent1.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                        parent1_copy.get_professor().set_schedule(new_day, new_time_slot)
                        remove_prof_sched = True

            if "c" in attribute:
                parent1_copy.set_course(item.get_course())
            if "p" in attribute:
                parent1_copy.set_professor(item.get_professor())

                new_time_slot = parent1_copy.get_class_slot_one(2)
                new_day = parent1_copy.get_class_slot_one(1)

                if not remove_prof_sched:
                    parent1.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                    parent1_copy.get_professor().set_schedule(new_day, new_time_slot)

                else:
                    curr_time_slot = parent1_copy.get_class_slot_one(2)
                    curr_day = parent1_copy.get_class_slot_one(1)

                    parent1.get_professor().get_schedule().remove_sched(curr_day, curr_time_slot)
                    parent1_copy.get_professor().set_schedule(curr_day, curr_time_slot)

            if "b" in attribute:
                parent1_copy.set_block(item.get_block())

                new_time_slot = parent1_copy.get_class_slot_one(2)
                new_day = parent1_copy.get_class_slot_one(1)

                if not remove_block_sched:
                    parent1.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                    parent1_copy.get_block().set_schedule(new_day, new_time_slot)
                else:
                    curr_time_slot = parent1_copy.get_class_slot_one(2)
                    curr_day = parent1_copy.get_class_slot_one(1)

                    parent1.get_block().get_schedule().remove_sched(curr_day, curr_time_slot)
                    parent1_copy.get_block().set_schedule(curr_day, curr_time_slot)

            gene_fitness = GeneFitness(parent1_copy)
            fitness_score = GeneticAlgorithm.fitness_value(gene_fitness)

            parent1.get_block().get_schedule().clear_and_add_new_sched(prev_day, orig_block_sched)
            parent1.get_professor().get_schedule().clear_and_add_new_sched(prev_day, orig_prof_sched)
            item.get_block().get_schedule().clear_and_add_new_sched(item_prev_day, item_orig_block_sched1)
            item.get_professor().get_schedule().clear_and_add_new_sched(item_prev_day, item_orig_prof_sched1)

            if different_day:
                parent1.get_block().get_schedule().clear_and_add_new_sched(item_prev_day, orig_block_sched2)
                parent1.get_professor().get_schedule().clear_and_add_new_sched(item_prev_day, orig_prof_sched2)
                item.get_block().get_schedule().clear_and_add_new_sched(prev_day, item_orig_block_sched2)
                item.get_professor().get_schedule().clear_and_add_new_sched(prev_day, item_orig_prof_sched2)

            if fitness_score == 1:
                # return parent2.append([item, 1])
                return item
            elif fitness_score > perfect_gene_score:
                print("found2 : ", item)
                perfect_gene_found = item
                perfect_gene_score = fitness_score

        # return the index of perfect gene for parent 1
        # return parent2.append([perfect_gene_found, perfect_gene_score])
        return perfect_gene_found

    @staticmethod
    def fitness_value(gene_fitness):
        if gene_fitness.first_year() is None:
            first_year = 1
        else:
            first_year = gene_fitness.first_year()

        fitness_score = (
                (
                        gene_fitness.classroom_capacity() +
                        gene_fitness.professor_work_load() +
                        gene_fitness.schedule_availability() +
                        gene_fitness.room_suitability() +
                        gene_fitness.course_slot_suitability() +
                        first_year +
                        gene_fitness.block_lunch() +
                        gene_fitness.professor_lunch() +
                        gene_fitness.maximum_working_hours() +
                        gene_fitness.prof_handled_course() +
                        gene_fitness.block_enrolled_course()
                ) / 12
        )

        return fitness_score

    def uniform_crossover(self, parent1, parent2):
        attribute = GeneticAlgorithm.check_attribute(parent1)

        parent1_copy = copy.copy(parent1)
        parent2_copy = copy.copy(parent2)

        prev_day = parent1_copy.get_class_slot_one(1)

        print("uniform crossover", )
        print("Parent 1")
        print("class slot : ", parent1.get_class_slot())
        print("Room sched: ", parent1.get_class_slot_one(0).get_schedule().get_schedule(prev_day))
        print("Prof : ", parent1.get_professor())
        print("Prof sched: ", parent1.get_professor().get_schedule().get_schedule(prev_day))
        print("Block : ", parent1.get_block())
        print("Block sched: ", parent1.get_block().get_schedule().get_schedule(prev_day))
        print("Course : ", parent1.get_course())

        remove_prof_sched = False
        remove_block_sched = False

        if 'cs' in attribute or 'r' in attribute or 'd' in attribute:
            print("cs uniform c")

            parent1.set_class_slot(parent2_copy.get_class_slot())
            parent2.set_class_slot(parent1_copy.get_class_slot())

            p1_prev_time_slot = parent1_copy.get_class_slot_one(2)
            p1_prev_day = parent1_copy.get_class_slot_one(1)
            p2_prev_time_slot = parent2_copy.get_class_slot_one(2)
            p2_prev_day = parent2_copy.get_class_slot_one(1)

            p1_new_time_slot = parent1.get_class_slot_one(2)
            p1_new_day = parent1.get_class_slot_one(1)
            p2_new_time_slot = parent2.get_class_slot_one(2)
            p2_new_day = parent2.get_class_slot_one(1)

            # prof sched will change
            # block sched will change
            parent1_copy.get_block().get_schedule().remove_sched(p1_prev_day, p1_prev_time_slot)
            parent1_copy.get_professor().get_schedule().remove_sched(p1_prev_day, p1_prev_time_slot)
            parent1.get_block().set_schedule(p1_new_day, p1_new_time_slot)
            parent1.get_professor().set_schedule(p1_new_day, p1_new_time_slot)

            parent2_copy.get_block().get_schedule().remove_sched(p2_prev_day, p2_prev_time_slot)
            parent2_copy.get_professor().get_schedule().remove_sched(p2_prev_day, p2_prev_time_slot)
            parent2.get_block().set_schedule(p2_new_day, p2_new_time_slot)
            parent2.get_professor().set_schedule(p2_new_day, p2_new_time_slot)

            print("prev p1 class slot : ", parent1_copy.get_class_slot())
            print("prev p1 Room sched: ", parent1_copy.get_class_slot_one(0).get_schedule().get_schedule(prev_day))
            print("prev p1 Prof sched: ", parent1_copy.get_professor().get_schedule().get_schedule(prev_day))
            print("prev p1 Block sched: ", parent1_copy.get_block().get_schedule().get_schedule(prev_day))

            print("prev p2 class slot : ", parent2_copy.get_class_slot())
            print("prev p2 Room sched: ", parent2_copy.get_class_slot_one(0).get_schedule().get_schedule(prev_day))
            print("prev p2 Prof sched: ", parent2_copy.get_professor().get_schedule().get_schedule(prev_day))
            print("prev p2 Block sched: ", parent2_copy.get_block().get_schedule().get_schedule(prev_day))

            print("new p1 class slot : ", parent1.get_class_slot())
            print("new p1 Room sched: ", parent1.get_class_slot_one(0).get_schedule().get_schedule(prev_day))
            print("new p1 Prof sched: ", parent1.get_professor().get_schedule().get_schedule(prev_day))
            print("new p1 Block sched: ", parent1.get_block().get_schedule().get_schedule(prev_day))

            print("new p2 class slot : ", parent2.get_class_slot())
            print("new p2 Room sched: ", parent2.get_class_slot_one(0).get_schedule().get_schedule(prev_day))
            print("new p2 Prof sched: ", parent2.get_professor().get_schedule().get_schedule(prev_day))
            print("new p2 Block sched: ", parent2.get_block().get_schedule().get_schedule(prev_day))

            remove_prof_sched = True
            remove_block_sched = True

        if 'p' in attribute:

            print("p uniform c")

            p1_prev_time_slot = parent1_copy.get_class_slot_one(2)
            p1_prev_day = parent1_copy.get_class_slot_one(1)
            p2_prev_time_slot = parent2_copy.get_class_slot_one(2)
            p2_prev_day = parent2_copy.get_class_slot_one(1)

            p1_new_time_slot = parent1.get_class_slot_one(2)
            p1_new_day = parent1.get_class_slot_one(1)
            p2_new_time_slot = parent2.get_class_slot_one(2)
            p2_new_day = parent2.get_class_slot_one(1)

            parent1.set_professor(parent2_copy.get_professor())
            parent2.set_professor(parent1_copy.get_professor())

            if not remove_prof_sched:
                parent1_copy.get_professor().get_schedule().remove_sched(p1_prev_day, p1_prev_time_slot)
                parent2_copy.get_professor().get_schedule().remove_sched(p2_prev_day, p2_prev_time_slot)

                parent1.get_professor().set_schedule(p1_new_day, p1_new_time_slot)
                parent2.get_professor().set_schedule(p2_new_day, p2_new_time_slot)
            else:
                p1_curr_time_slot = parent1.get_class_slot_one(2)
                p1_curr_day = parent1.get_class_slot_one(1)
                p2_curr_time_slot = parent2.get_class_slot_one(2)
                p2_curr_day = parent2.get_class_slot_one(1)

                parent1_copy.get_professor().get_schedule().remove_sched(p1_curr_day, p1_curr_time_slot)
                parent2_copy.get_professor().get_schedule().remove_sched(p2_curr_day, p2_curr_time_slot)
                parent1.get_professor().set_schedule(p1_curr_day, p1_curr_time_slot)
                parent2.get_professor().set_schedule(p2_curr_day, p2_curr_time_slot)

        if 'b' in attribute:

            print("b uniform c")

            p1_prev_time_slot = parent1_copy.get_class_slot_one(2)
            p1_prev_day = parent1_copy.get_class_slot_one(1)
            p2_prev_time_slot = parent2_copy.get_class_slot_one(2)
            p2_prev_day = parent2_copy.get_class_slot_one(1)

            p1_new_time_slot = parent1.get_class_slot_one(2)
            p1_new_day = parent1.get_class_slot_one(1)
            p2_new_time_slot = parent2.get_class_slot_one(2)
            p2_new_day = parent2.get_class_slot_one(1)

            parent1.set_block(parent2_copy.get_block())
            parent2.set_block(parent1_copy.get_block())

            if not remove_block_sched:
                parent1_copy.get_block().get_schedule().remove_sched(p1_prev_day, p1_prev_time_slot)
                parent2_copy.get_block().get_schedule().remove_sched(p2_prev_day, p2_prev_time_slot)

                parent1.get_block().set_schedule(p1_new_day, p1_new_time_slot)
                parent2.get_block().set_schedule(p2_new_day, p2_new_time_slot)
            else:
                p1_curr_time_slot = parent1.get_class_slot_one(2)
                p1_curr_day = parent1.get_class_slot_one(1)
                p2_curr_time_slot = parent2.get_class_slot_one(2)
                p2_curr_day = parent2.get_class_slot_one(1)

                parent1_copy.get_block().get_schedule().remove_sched(p1_curr_day, p1_curr_time_slot)
                parent2_copy.get_block().get_schedule().remove_sched(p2_curr_day, p2_curr_time_slot)
                parent1.get_block().set_schedule(p1_curr_day, p1_curr_time_slot)
                parent2.get_block().set_schedule(p2_curr_day, p2_curr_time_slot)

        if 'c' in attribute:
            print("c uniform c")

            parent1.set_course(parent2_copy.get_course())
            parent2.set_course(parent1_copy.get_course())

        parent1_fitness = GeneFitness(parent1)
        parent1_fitness_score = GeneticAlgorithm.fitness_value(parent1_fitness)
        parent1.add_fitness_score(parent1_fitness_score)

        parent2_fitness = GeneFitness(parent2)
        parent2_fitness_score = GeneticAlgorithm.fitness_value(parent2_fitness)
        parent2.add_fitness_score(parent2_fitness_score)

        print("p1 before crossover : ", parent1_copy.get_fitness_score())
        print("p2 before crossover : ", parent2_copy.get_fitness_score())
        print("p1 after crossover : ", parent1.get_fitness_score())
        print("p2 after crossover : ", parent2.get_fitness_score(), "index : ",
              self.curr_chromosome.get_position(parent2))

        if parent1.get_fitness_score() == 1 and parent2.get_fitness_score() == 1:
            print("luh")
            parent2 = self.curr_chromosome.get_genes().pop(self.curr_chromosome.get_position(parent2))
            self.add_to_new_chromosome(parent1)
            self.add_to_new_chromosome(parent2)
            self.chromosome_size = self.chromosome_size - 1
            print("update size 11 : ", self.get_chromosome_size())
        elif parent1.get_fitness_score() == 1 and parent2.get_fitness_score() != 1:
            self.add_to_new_chromosome(parent1)
        elif parent1.get_fitness_score() != 1 and parent2.get_fitness_score() == 1:
            print("luh")
            parent2 = self.curr_chromosome.get_genes().pop(self.curr_chromosome.get_position(parent2))
            self.add_to_new_chromosome(parent2)
            self.curr_chromosome.add_gene(parent1)
            print("update size 22 : ", self.get_chromosome_size())
        elif parent1.get_fitness_score() > parent2.get_fitness_score():
            print("leh")
            child = self.mutation(parent1)
            self.add_to_new_chromosome(child)
        elif parent1.get_fitness_score() < parent2.get_fitness_score():
            print("lih")
            parent2_copy = copy.copy(parent2)
            self.curr_chromosome.pop_gene(self.curr_chromosome.get_position(parent2))
            self.curr_chromosome.add_gene(parent1)
            child = self.mutation(parent2_copy)
            self.add_to_new_chromosome(child)
        else:
            child = self.mutation(parent1)
            self.add_to_new_chromosome(child)

        print("NEW CHROMOSOME LEN : ", len(self.new_chromosome.get_genes()))

    def mutation(self, child):

        remove_prof_sched = False
        remove_block_sched = False
        random_gene_picked = False
        child_copy = copy.copy(child)
        attribute = self.check_attribute(child)

        prev_time_slot = child_copy.get_class_slot_one(2)
        prev_day = child_copy.get_class_slot_one(1)

        if self.chromosome_size != 0:

            random_num = np.random.randint(self.chromosome_size)
            random_gene = self.curr_chromosome.get_specific_gene(random_num)
            rg_copy = copy.copy(random_gene)

            print("mutatiooon")
            print("child id : ", child)

            prev_time_slot_rg = rg_copy.get_class_slot_one(2)
            prev_day_rg = rg_copy.get_class_slot_one(1)

            print("prev cs rg", random_gene.get_class_slot())

            print("RANDOM GENE ORIGINAL")
            print("rg block sched : ", random_gene.get_block().get_schedule().get_schedule(prev_day_rg))
            print("rg prof sched : ", random_gene.get_block().get_schedule().get_schedule(prev_day_rg))

            if 'cs' in attribute or 'r' in attribute or 'd' in attribute:

                print("mutation cs")

                if len(self.class_slot_population) > 0:

                    random_class_slot = self.class_slot_population.pop(
                        random.randrange(len(self.class_slot_population)))

                    print("class slot size : ", len(self.class_slot_population))

                    child_class_slot = child.get_class_slot()

                    self.class_slot_population.append(child_class_slot)

                    child.set_class_slot(random_class_slot)

                    print("class slot size1 : ", len(self.class_slot_population))
                else:
                    child_class_slot = child.get_class_slot()
                    child.set_class_slot(random_gene.get_class_slot())
                    random_gene.set_class_slot(child_class_slot)
                    random_gene_picked = True

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)

                if random_gene_picked:
                    print("rg picked")
                    rg_new_time_slot = random_gene.get_class_slot_one(2)
                    rg_new_day = random_gene.get_class_slot_one(1)
                    rg_copy.get_block().get_schedule().remove_sched(prev_day_rg, prev_time_slot_rg)
                    rg_copy.get_professor().get_schedule().remove_sched(prev_day_rg, prev_time_slot_rg)
                    random_gene.get_block().set_schedule(rg_new_day, rg_new_time_slot)
                    random_gene.get_professor().set_schedule(rg_new_day, rg_new_time_slot)

                child_copy.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                child_copy.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                child.get_block().set_schedule(new_day, new_time_slot)
                child.get_professor().set_schedule(new_day, new_time_slot)

                print("cs mutation")
                print("prev time slot: ", prev_time_slot)
                print("prev block sched : ", child_copy.get_block().get_schedule().get_schedule(prev_day))
                print("prev prof sched : ", child_copy.get_professor().get_schedule().get_schedule(prev_day))
                print("new time slot: ", new_time_slot)
                print("new block sched : ", child.get_block().get_schedule().get_schedule(new_day))
                print("new prof sched : ", child.get_professor().get_schedule().get_schedule(new_day))

                remove_block_sched = True
                remove_prof_sched = True

            if 'p' in attribute:
                print("mutation p")

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)
                rg_new_time_slot = random_gene.get_class_slot_one(2)
                rg_new_day = random_gene.get_class_slot_one(1)

                child_prof = child.get_professor()
                print("mutation prof : ", child_prof)

                child.set_professor(random_gene.get_professor())
                random_gene.set_professor(child_prof)

                if not remove_prof_sched:
                    child_copy.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                    rg_copy.get_professor().get_schedule().remove_sched(prev_day_rg, prev_time_slot_rg)

                    child.get_professor().set_schedule(new_day, new_time_slot)
                    random_gene.get_professor().set_schedule(rg_new_day, rg_new_time_slot)

                else:
                    curr_time_slot = child.get_class_slot_one(2)
                    curr_day = child.get_class_slot_one(1)
                    rg_curr_time_slot = random_gene.get_class_slot_one(2)
                    rg_curr_day = random_gene.get_class_slot_one(1)

                    if random_gene_picked:
                        child_copy.get_professor().get_schedule().remove_sched(curr_day, curr_time_slot)
                        rg_copy.get_professor().get_schedule().remove_sched(rg_curr_day, rg_curr_time_slot)
                        child.get_professor().set_schedule(curr_day, curr_time_slot)
                        random_gene.get_professor().set_schedule(rg_curr_day, rg_curr_time_slot)
                    else:
                        child_copy.get_professor().get_schedule().remove_sched(curr_day, curr_time_slot)
                        rg_copy.get_professor().get_schedule().remove_sched(prev_day_rg, prev_time_slot_rg)
                        child.get_professor().set_schedule(curr_day, curr_time_slot)
                        random_gene.get_professor().set_schedule(rg_curr_day, rg_curr_time_slot)

            if 'b' in attribute:
                print("mutation b")

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)
                rg_new_time_slot = random_gene.get_class_slot_one(2)
                rg_new_day = random_gene.get_class_slot_one(1)

                child_block = child.get_block()
                random_gene_block = random_gene.get_block()

                child.set_block(random_gene_block)
                random_gene.set_block(child_block)

                print("RANDOM GENE BLOCK MUTATION")
                print("rg block sched : ", rg_copy.get_block().get_schedule().get_schedule(prev_day_rg))
                print("rg block sched1 : ", random_gene.get_block().get_schedule().get_schedule(rg_new_day))

                if not remove_block_sched:
                    child_copy.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                    rg_copy.get_block().get_schedule().remove_sched(prev_day_rg, prev_time_slot_rg)

                    child.get_block().set_schedule(new_day, new_time_slot)
                    random_gene.get_block().set_schedule(rg_new_day, rg_new_time_slot)
                else:

                    curr_time_slot = child.get_class_slot_one(2)
                    curr_day = child.get_class_slot_one(1)
                    rg_curr_time_slot = random_gene.get_class_slot_one(2)
                    rg_curr_day = random_gene.get_class_slot_one(1)

                    if random_gene_picked:
                        child_copy.get_block().get_schedule().remove_sched(curr_day, curr_time_slot)
                        rg_copy.get_block().get_schedule().remove_sched(rg_curr_day, rg_curr_time_slot)
                        child.get_block().set_schedule(curr_day, curr_time_slot)
                        random_gene.get_block().set_schedule(rg_curr_day, rg_curr_time_slot)
                    else:
                        child_copy.get_block().get_schedule().remove_sched(curr_day, curr_time_slot)
                        rg_copy.get_block().get_schedule().remove_sched(prev_day_rg, prev_time_slot_rg)
                        child.get_block().set_schedule(curr_day, curr_time_slot)
                        random_gene.get_block().set_schedule(rg_curr_day, rg_curr_time_slot)

            if 'c' in attribute:
                child_course = child.get_course()
                print("mutation course : ", child_course)
                child.set_course(random_gene.get_course())
                random_gene.set_course(child_course)
        else:
            if 'cs' in attribute or 'r' in attribute or 'd' in attribute:

                print("mutation cs")

                if len(self.class_slot_population) > 0:

                    random_class_slot = self.class_slot_population.pop(
                        random.randrange(len(self.class_slot_population)))

                    child_class_slot = child.get_class_slot()
                    self.class_slot_population.append(child_class_slot)
                    child.set_class_slot(random_class_slot)

                    new_time_slot = child.get_class_slot_one(2)
                    new_day = child.get_class_slot_one(1)

                    child_copy.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                    child_copy.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                    child.get_block().set_schedule(new_day, new_time_slot)
                    child.get_professor().set_schedule(new_day, new_time_slot)

        return child

    def add_to_new_chromosome(self, gene):
        self.new_chromosome.add_gene(gene)

    def get_chromosome_size(self):
        return len(self.curr_chromosome.get_genes())

    def add_new_chromosome(self):
        new_chromosome = Chromosome()

        GeneticAlgoPopulation.chromosome.append(self.new_chromosome)
        self.new_chromosome = new_chromosome
        self.chromosome_size = 0


class Distributor:
    def __init__(self, items, num_threads, attribute, parent1):
        self.items = items
        self.num_threads = num_threads
        self.parent1 = parent1
        self.parent2 = []
        self.attribute = attribute
        self.thread_items = [[] for _ in range(num_threads)]
        self.distribute_items()

    def distribute_items(self):
        num_items = len(self.items)
        items_per_thread = num_items // self.num_threads
        remainder = num_items % self.num_threads

        index = 0
        for i in range(self.num_threads):
            if remainder > 0:
                num_assigned = items_per_thread + 1
                remainder -= 1
            else:
                num_assigned = items_per_thread

            self.thread_items[i] = self.items[index:index + num_assigned]
            index += num_assigned

    def start_threads(self):
        threads = []
        for i in range(self.num_threads):
            thread = threading.Thread(target=GeneticAlgorithm.parent_two, args=(
                self.thread_items[i], self.attribute, self.parent1, self.parent2))
            thread.start()
            threads.append(thread)

        for thread in threads:
            thread.join()

        return self.parent2


class Chromosome:
    def __init__(self):
        self.genes = []
        self.fitness_value = 0

    def add_gene(self, gene):
        self.genes.append(gene)

    def get_genes(self):
        return self.genes

    def get_specific_gene(self, index):
        return self.genes[index]

    def pop_gene(self, index):
        self.genes.pop(index)

    def get_position(self, gene):
        return self.genes.index(gene)

    def set_fitness_value(self, value):
        self.fitness_value = value

    def get_fitness_value(self):
        return self.fitness_value

    def clear_chromosome(self):
        self.genes.clear()


class Gene:
    def __init__(self, class_slot, course, block, professor):
        self.class_slot = class_slot
        self.course = course
        self.block = block
        self.professor = professor
        self.fitness_score = 0

        self.classroom_capacity_fitness = 0
        self.professor_work_load_fitness = 0
        self.room_availability_fitness = 0
        self.block_availability_fitness = 0
        self.room_suitability_fitness = 0
        self.course_slot_suitability_fitness = 0
        self.first_year_fitness = 0
        self.prof_lunch_break_fitness = 0
        self.block_lunch_break_fitness = 0
        self.working_hours_fitness = 0
        self.prof_handled_course_fitness = 0
        self.block_enrolled_course_fitness = 0

    def get_class_slot(self):
        return self.class_slot

    def get_class_slot_one(self, index):
        return self.class_slot[index]

    def get_class_slot_two(self, index1, index2):
        return self.class_slot[index1][index2]

    def get_course(self):
        return self.course

    def get_block(self):
        return self.block

    def get_professor(self):
        return self.professor

    def set_class_slot(self, value):
        self.class_slot = value

    def set_course(self, value):
        self.course = value

    def set_block(self, value):
        self.block = value

    def set_professor(self, value):
        self.professor = value

    def add_fitness_score(self, fitness):
        self.fitness_score = fitness

    def get_fitness_score(self):
        return self.fitness_score

    def get_classroom_capacity_fitness(self):
        return self.classroom_capacity_fitness

    def get_professor_work_load_fitness(self):
        return self.professor_work_load_fitness

    def get_room_availability_fitness(self):
        return self.room_availability_fitness

    def get_block_availability_fitness(self):
        return self.block_availability_fitness

    def get_room_suitability_fitness(self):
        return self.room_suitability_fitness

    def get_first_year_fitness(self):
        return self.first_year_fitness

    def get_prof_lunch_break_fitness(self):
        return self.prof_lunch_break_fitness

    def get_block_lunch_break_fitness(self):
        return self.block_lunch_break_fitness

    def get_working_hours_fitness(self):
        return self.working_hours_fitness

    def get_course_slot_suitability_fitness(self):
        return self.course_slot_suitability_fitness

    def get_prof_handled_course_fitness(self):
        return self.prof_handled_course_fitness

    def get_block_enrolled_course_fitness(self):
        return self.block_enrolled_course_fitness

    def set_classroom_capacity_fitness(self, value):
        self.classroom_capacity_fitness = value

    def set_professor_work_load_fitness(self, value):
        self.professor_work_load_fitness = value

    def set_room_availability_fitness(self, value):
        self.room_availability_fitness = value

    def set_block_availability_fitness(self, value):
        self.block_availability_fitness = value

    def set_room_suitability_fitness(self, value):
        self.room_suitability_fitness = value

    def set_first_year_fitness(self, value):
        self.first_year_fitness = value

    def set_prof_lunch_break_fitness(self, value):
        self.prof_lunch_break_fitness = value

    def set_block_lunch_break_fitness(self, value):
        self.block_lunch_break_fitness = value

    def set_working_hours_fitness(self, value):
        self.working_hours_fitness = value

    def set_course_slot_suitability_fitness(self, value):
        self.course_slot_suitability_fitness = value

    def set_prof_handled_course_fitness(self, value):
        self.prof_handled_course_fitness = value

    def set_block_enrolled_course_fitness(self, value):
        self.block_enrolled_course_fitness = value


class GeneFitness:

    def __init__(self, gene):
        self.gene = gene

    def classroom_capacity(self):
        if self.gene.block.get_enrolled_students() <= self.gene.get_class_slot_one(0).get_capacity():
            self.gene.set_classroom_capacity_fitness(1)
            return 1
        else:
            self.gene.set_classroom_capacity_fitness(0)
            return 0

    def professor_work_load(self):
        day = self.gene.get_class_slot_one(1)
        schedule = self.gene.professor.get_schedule().get_schedule(day)

        if MEMOIZATION:
            for sublist in PWL_CACHE:
                if sublist[0] == [self.gene.get_professor(), day]:
                    print("PWL MEMO0 : ", sublist[0])
                    print("PWL MEMO : ", sublist[1])
                    self.gene.set_professor_work_load_fitness(sublist[1])
                    print("PWL FIT : ", self.gene.get_professor_work_load_fitness())
                    return sublist[1]
        print("PWL RUNNING")
        result = has_overlap_and_multiple_copies(self.gene.get_class_slot_one(2), schedule)

        if result:
            self.gene.set_professor_work_load_fitness(0)
            if MEMOIZATION:
                PWL_CACHE.append([[self.gene.get_professor(), day], 0])
            return 0
        else:
            self.gene.set_professor_work_load_fitness(1)
            if MEMOIZATION:
                PWL_CACHE.append([[self.gene.get_professor(), day], 1])
            return 1

    def schedule_availability(self):
        day = self.gene.get_class_slot_one(1)
        # print("dayday: ", day)
        room_schedule = self.gene.class_slot[0].get_schedule().get_schedule(day)
        block_schedule = self.gene.block.get_schedule().get_schedule(day)

        room_availability = has_overlap_and_multiple_copies(self.gene.get_class_slot_one(2), room_schedule)
        block_availability = has_overlap_and_multiple_copies(self.gene.get_class_slot_one(2), block_schedule)

        if room_availability:
            self.gene.set_room_availability_fitness(0)
        else:
            self.gene.set_room_availability_fitness(1)

        if block_availability:
            self.gene.set_block_availability_fitness(0)
        else:
            self.gene.set_block_availability_fitness(1)

        if not room_availability and not block_availability:
            return 2
        elif not room_availability or not block_availability:
            return 1
        elif room_availability and block_availability:
            return 0

    def room_suitability(self):
        room = self.gene.get_class_slot_one(0)
        course = self.gene.course
        available_rooms = course.get_available_rooms()

        # print("course.get_available_rooms() : ", course.get_available_rooms())

        if room.get_code() in available_rooms:
            self.gene.set_room_suitability_fitness(1)
            return 1
        else:
            self.gene.set_room_suitability_fitness(0)
            return 0

    def course_slot_suitability(self):
        slot = self.gene.get_class_slot_one(2)

        # print("length : ", len(slot) - 1, " hour : ", self.gene.course.get_hour())
        if len(slot) - 1 == self.gene.course.get_hour():
            self.gene.set_course_slot_suitability_fitness(1)
            return 1
        else:
            self.gene.set_course_slot_suitability_fitness(0)
            return 0

    def first_year(self):

        block = Data.find_block(self.gene.block)

        # print("year : ", block.get_year())
        if block.get_year() == 1:
            if self.gene.get_class_slot_one(1) != 'S':
                self.gene.set_first_year_fitness(1)
                return 1
            else:
                self.gene.set_first_year_fitness(0)
                return 0
        else:
            self.gene.set_first_year_fitness(2)
            return None

    def professor_lunch(self):
        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        if MEMOIZATION:
            for sublist in PL_CACHE:
                if sublist[0] == [self.gene.get_professor(), day]:
                    PL_CACHE.append([[self.gene.get_professor(), day], sublist[1]])
                    return sublist[1]

        # if true then no lunch
        if has_no_lunch_break(professor_schedule):
            self.gene.set_prof_lunch_break_fitness(0)
            if MEMOIZATION:
                PL_CACHE.append([[self.gene.get_professor(), day], 0])
            return 0
        else:
            self.gene.set_prof_lunch_break_fitness(1)
            if MEMOIZATION:
                PL_CACHE.append([[self.gene.get_professor(), day], 1])
            return 1

    def block_lunch(self):
        block_sched_day = self.gene.get_class_slot_one(1)
        block_schedule = self.gene.block.schedule.get_schedule(block_sched_day)

        if MEMOIZATION:
            for sublist in BL_CACHE:
                if sublist[0] == [self.gene.get_block(), day]:
                    BL_CACHE.append([[self.gene.get_block(), day], sublist[1]])
                    return sublist[1]

        # if true then no lunch
        if has_no_lunch_break(block_schedule):
            self.gene.set_block_lunch_break_fitness(0)
            if MEMOIZATION:
                BL_CACHE.append([[self.gene.get_block(), day], 0])
            return 0
        else:
            self.gene.set_block_lunch_break_fitness(1)
            if MEMOIZATION:
                BL_CACHE.append([[self.gene.get_block(), day], 1])
            return 1

    def maximum_working_hours(self):
        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        if MEMOIZATION:
            for sublist in WH_CACHE:
                if sublist[0] == [self.gene.get_professor(), day]:
                    if sublist[1] != 1:
                        self.gene.set_working_hours_fitness(0)
                    else:
                        self.gene.set_working_hours_fitness(1)
                    return sublist[1]

        total = 0
        for sublist in professor_schedule:
            length = len(sublist) - 1
            total += length

        if total > 8:
            total = 8 - total
            self.gene.set_working_hours_fitness(0)
            if MEMOIZATION:
                WH_CACHE.append([[self.gene.get_professor(), day, total]])
        else:
            self.gene.set_working_hours_fitness(1)
            total = 1
            if MEMOIZATION:
                WH_CACHE.append([[self.gene.get_professor(), day], 1])

        return total

    def prof_handled_course(self):
        course_code = self.gene.course.get_code()

        if course_code in self.gene.professor.get_course_code_handle():
            self.gene.set_prof_handled_course_fitness(1)
            return 1
        else:
            self.gene.set_prof_handled_course_fitness(0)
            return 0

    def block_enrolled_course(self):
        course_code = self.gene.course.get_code()

        if course_code in self.gene.block.get_courses():
            self.gene.set_block_enrolled_course_fitness(1)
            return 1
        else:
            self.gene.set_block_enrolled_course_fitness(0)
            return 0


class ClassRoom:

    def __init__(self, code, capacity, type_of_room):
        self.code = code
        self.capacity = capacity
        self.type_of_room = type_of_room
        self.schedule = Schedule()

    def get_code(self):
        return self.code

    def get_capacity(self):
        return self.capacity

    def get_schedule(self):
        return self.schedule

    def get_type_of_room(self):
        return self.type_of_room

    def set_schedule(self, day, time_slot):
        self.schedule.set_schedule(day, time_slot)


class LabCourse:

    def __init__(self, code, hour, available_rooms):
        self.code = code
        self.hour = hour
        self.available_rooms = available_rooms

    def get_code(self):
        return self.code

    def get_hour(self):
        return self.hour

    def get_available_rooms(self):
        return self.available_rooms

    def get_year(self):
        return self.year


class LectureCourse:

    def __init__(self, code, hour, available_rooms):
        self.code = code
        self.hour = hour
        self.available_rooms = available_rooms

    def get_code(self):
        return self.code

    def get_hour(self):
        return self.hour

    def get_available_rooms(self):
        return self.available_rooms

    def get_year(self):
        return self.year


class Professor:

    def __init__(self, prof_id, course_code_handle):
        self.prof_id = prof_id
        self.course_code_handle = course_code_handle
        self.schedule = Schedule()

    def get_prof_id(self):
        return self.prof_id

    def get_course_code_handle(self):
        return self.course_code_handle

    def get_schedule(self):
        return self.schedule

    def set_schedule(self, day, time_slot):
        self.schedule.set_schedule(day, time_slot)


class Block:

    def __init__(self, block_code, enrolled_students, courses, year):
        self.block_code = block_code
        self.enrolled_students = enrolled_students
        self.courses = courses
        self.year = year
        self.schedule = Schedule()

    def get_block_code(self):
        return self.block_code

    def get_enrolled_students(self):
        return self.enrolled_students

    def get_courses(self):
        return self.courses

    def get_schedule(self):
        return self.schedule

    def set_schedule(self, day, time_slot):
        self.schedule.set_schedule(day, time_slot)

    def get_year(self):
        return self.year


class Population:
    def __init__(self):
        self.population = []
        self.courses = []
        self.blocks = []
        self.professors = []
        self.class_slots = []

    def initialize_population(self, meetingTime, meetingDay):
        self.generate_course_block_population()
        self.generate_prof_population()
        self.generate_class_slot_population(meetingTime, meetingDay)

    def generate_course_block_population(self):
        data = Data()
        blocks = data.block

        for block in blocks:
            for course in block.get_courses():
                for data_lec_course in Data.lec_course:
                    if course == data_lec_course.get_code():
                        self.courses.append(Data.find_lec_course(course))
                        self.blocks.append(block)
                for data_lab_course in Data.lab_course:
                    if course == data_lab_course.get_code():
                        self.courses.append(Data.find_lab_course(course))
                        self.blocks.append(block)

    def generate_prof_population(self):
        data = Data()

        for prof in data.professor:
            for course in prof.get_course_code_handle():
                for data_lec_course in Data.lec_course:
                    if course == data_lec_course.get_code():
                        self.professors.append(prof)
                for data_lab_course in Data.lab_course:
                    if course == data_lab_course.get_code():
                        self.professors.append(prof)

    def generate_class_slot_population(self, meetingTime, day):
        self.class_slots.extend(labRoomMapping(meetingTime, day))
        self.class_slots.extend(regRoomMapping(day))

    def get_population(self):
        return self.population

    def get_courses(self):
        return self.courses

    def get_blocks(self):
        return self.blocks

    def get_class_slots(self):
        return self.class_slots

    def get_professors(self):
        return self.professors

    def add_class_slot(self, class_slot):
        self.class_slots.append(class_slot)


def labRoomMapping(meetingTime: List[str], meetingDays: List[str]) -> List[List[str]]:
    """
        This function maps the available laboratory rooms with the given meeting time and meeting days.

        Arguments:
        meetingTime -- A list of meeting hours in the format "HH".
        meetingDays -- A list of meeting days.

        Returns:
        lab_map -- A list of lists, where each inner list represents the mapping of a laboratory room with a day and three hours of meeting time.

        Time complexity:  O(len(labRoomList) * (len(meetingTime)/3) * len(meetingDays)) or O(n^3)
        This is because, for each laboratory room in labRoomList, it loops through each day in meetingDays and each three-hour time slot in meetingTime. The length of each time slot is determined by the range function, where it steps by 3. The time complexity is proportional to the number of iterations required to map each laboratory room with its corresponding day and time slot.
        """

    labRoom = Data()
    lab_map = [[lab, day, meetingTime[i:i + 4]]
               for lab in labRoom.lab_room
               for day in meetingDays
               for i in range(0, len(meetingTime), 3)
               if len(meetingTime[i:i + 4]) == 4]

    return lab_map


def regRoomMapping(meetingDays: List[str]) -> List[List[str]]:
    reg_room = Data()
    pattern = [["07", "08", "09"], ["07", "08", "09", "10"], ["09", "10", "11"], ["09", "10", "11", "12"],
               ["10", "11", "12"], ["10", "11", "12", "13"],
               ["11", "12", "13"], ["11", "12", "13", "14"], ["12", "13", "14"], ["12", "13", "14", "15"],
               ["13", "14", "15"], ["13", "14", "15", "16"],
               ["14", "15", "16"], ["14", "15", "16", "17"], ["15", "16", "17"], ["16", "17", "18", "19"],
               ["17", "18", "19"]]

    reg_room_map = [[room, day, meetingTime]
                    for room in reg_room.reg_room
                    for day in meetingDays
                    for meetingTime in pattern
                    ]

    return reg_room_map


def has_overlap_and_multiple_copies(time_slots, pattern):
    # Check for more than 1 copy of pattern in time slots
    count = 0
    for p in pattern:
        if p == time_slots:
            count += 1
    if count > 1:
        return True

    # Remove duplicates of pattern
    pattern_copy = list(filter(lambda x: x != time_slots, pattern))

    if len(pattern_copy) > 0:
        # check for overlap in pattern
        if len(time_slots) > 3:
            for p in pattern_copy:
                if set(time_slots[1:4]).issubset(set(p)):
                    return True
                if set(time_slots[0:4]).issubset(set(p)):
                    return True
                if set(time_slots[0:3]).issubset(set(p)):
                    return True
                if set(time_slots[0:2]).issubset(set(p)):
                    return True
                elif set(time_slots[2:4]).issubset(set(p)):
                    return True
        else:
            for p in pattern_copy:
                if set(time_slots[0:3]).issubset(set(p)):
                    return True
                if set(time_slots[0:2]).issubset(set(p)):
                    return True
                if set(time_slots[1:3]).issubset(set(p)):
                    return True

    return False


def available_hours_in_schedule(schedule):
    merged = []
    for slot in schedule:
        # Convert slot to a set of integers
        slot_set = set(map(int, slot))
        # Check if there are any merged slots that overlap with this slot
        overlap_found = False
        for merged_slot in merged:
            # Convert merged_slot to a set of integers
            merged_slot_set = set(map(int, merged_slot))
            # Check if there is any overlap between slot_set and merged_slot_set
            if len(slot_set & merged_slot_set) > 0:
                # Merge the two sets and update merged_slot
                merged_slot |= slot_set
                overlap_found = True
                break
        # If no overlap was found, add slot to merged list as a new set
        if not overlap_found:
            merged.append(slot_set)
    # Convert the sets in merged back to lists of strings
    m = [list(map(str, s)) for s in merged]

    return 12 - (sum(map(len, m)) - len(m))


def has_no_lunch_break(pb_schedule):
    # Get all the sublists that contain a number between 10 and 14
    schedules = [schedule for schedule in pb_schedule if any(10 <= int(time) <= 14 for time in schedule)]

    # Create a new list with a sublist that contains only two elements
    relevant_sublists = []
    for schedule in schedules:
        for i in range(len(schedule) - 1):
            sublist = [int(schedule[i]), int(schedule[i + 1])]
            if 10 <= sublist[0] <= 14 and 10 <= sublist[1] <= 14:
                relevant_sublists.append(sublist)

    # Sort the array
    sorted_sublists = sorted(relevant_sublists)

    return set([tuple(sublist) for sublist in [[10, 11], [11, 12], [12, 13], [13, 14]]]).issubset(
        [tuple(sublist) for sublist in sorted_sublists])


if __name__ == '__main__':
    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F", "S"]

    Data.add_room(ClassRoom("SL1", 50, "Laboratory"))
    Data.add_room(ClassRoom("SL2", 50, "Laboratory"))
    Data.add_room(ClassRoom("SL3", 50, "Laboratory"))
    Data.add_room(ClassRoom("CISCO", 50, "Laboratory"))
    Data.add_room(ClassRoom("101", 50, "Regular"))
    Data.add_room(ClassRoom("102", 50, "Regular"))
    Data.add_room(ClassRoom("103", 50, "Regular"))
    Data.add_room(ClassRoom("104", 50, "Regular"))

    # year = 1 if first year else 0
    Data.add_lab_course(LabCourse("CS 121", 3, ["Sl1", "SL2", "SL3"]))
    Data.add_lab_course(LabCourse("CS 222", 3, ["Sl1", "SL2", "SL3", "CISCO"]))
    Data.add_lab_course(LabCourse("CS 221", 3, ["Sl1", "SL2"]))
    Data.add_lab_course(LabCourse("CS 322", 3, ["Sl1", "SL2", "SL3"]))

    Data.add_lab_course(LabCourse("CS 121", 2, ["Sl1", "SL2", "SL3"]))
    Data.add_lab_course(LabCourse("CS 222", 2, ["Sl1", "SL2", "SL3", "CISCO"]))
    Data.add_lab_course(LabCourse("CS 221", 2, ["Sl1", "SL2"]))
    Data.add_lab_course(LabCourse("CS 322", 2, ["Sl1", "SL2", "SL3"]))

    Data.add_lec_course(LectureCourse("MATH 111", 3, ["101", "102", "103", "104"]))
    Data.add_lec_course(LectureCourse("GEd 105", 3, ["101", "102", "103", "104"]))
    Data.add_lec_course(LectureCourse("GEd 107", 3, ["101", "102", "103", "104"]))
    Data.add_lec_course(LectureCourse("CS 323", 3, ["101", "102", "103", "104"]))

    Data.add_block(Block("101", 45, ["CS 121", "CS 222", "MATH 111"], 1))
    Data.add_block(Block("102", 45, ["CS 121", "CS 221", "MATH 111"], 0))
    Data.add_block(Block("201", 45, ["CS 322", "CS 323", "GEd 107"], 0))
    Data.add_block(Block("202", 45, ["GEd 105", "CS 222", "CS 121"], 0))
    Data.add_block(Block("203", 45, ["GEd 105", "CS 222", "CS 121"], 0))
    Data.add_block(Block("301", 45, ["CS 322", "CS 323", "GEd 107"], 0))

    Data.add_prof(Professor("P1", ["CS 121", "CS 222", "MATH 111"]))
    Data.add_prof(Professor("P2", ["CS 121", "CS 221", "CS 323"]))
    Data.add_prof(Professor("P3", ["MATH 111", "CS 322", "GEd 107"]))
    Data.add_prof(Professor("P4", ["GEd 105", "CS 322", "CS 222", "CS 323"]))
    Data.add_prof(Professor("P5", ["CS 121", "GEd 107"]))
    Data.add_prof(Professor("P6", ["CS 121", "CS 222", "GEd 105"]))

    population = Population()
    population.initialize_population(meetingTime, day)

    genetic_algo = GeneticAlgorithm()
    genetic_algo.encode_chromosome(population)

    not_perfect_schedule = True

    while not_perfect_schedule:

        chromosome_fit = genetic_algo.calculate_fit()
        new_chromosome = True
        MEMOIZATION = False

        if chromosome_fit == 1:
            not_perfect_schedule = False
        else:
            while new_chromosome:
                parent1 = genetic_algo.roulette_wheel()

                if parent1 is not None:

                    parent2 = genetic_algo.hill_climbing(parent1)

                    if parent2 is None:
                        genetic_algo.add_to_new_chromosome(parent1)
                    elif parent2 == 0:
                        genetic_algo.mutation(parent1)
                    else:
                        genetic_algo.uniform_crossover(parent1, parent2)
                else:
                    print("PARENT 1 IS NONE")

                print("CHROMOSOME SIZE : ", genetic_algo.get_chromosome_size())

                if genetic_algo.get_chromosome_size() == 0:
                    print("CHROMOSOME IS EMPTY")
                    genetic_algo.add_new_chromosome()
                    WH_CACHE.clear()
                    PL_CACHE.clear()
                    BL_CACHE.clear()
                    PWL_CACHE.clear()
                    MEMOIZATION = True
                    new_chromosome = False
