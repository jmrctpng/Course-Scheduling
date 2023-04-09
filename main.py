from tabulate import tabulate
from typing import List
import random
import time
import numpy as np
import copy

PROF_MUTATE = None
PROF_MUTATE1 = None
BLOCK_MUTATE = None
BLOCK_MUTATE1 = None
WH_CACHE = []
PL_CACHE = []
BL_CACHE = []
PWL_CACHE = []
MEMOIZATION = True

LAB_COURSE = ['CS 121', 'CS 111', 'CS 131']
REG_COURSE = ['FILI 102', 'GEd 105', 'GEd 108', 'MATH 401', 'MATH 407', 'MATH 111', 'GEd 109']
PROFESSOR = []
BLOCK = []


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
    def find_lec_course(course_id):

        for course in Data.lec_course:
            if course.get_code() == course_id:
                return course

        return None

    @staticmethod
    def find_lab_course(course_id):

        course1 = None
        course2 = None

        for course in Data.lab_course:
            if course.get_code() == course_id:

                if course1 is None:
                    course1 = course
                else:
                    course2 = course

                if course1 is not None and course2 is not None:
                    return course1, course2

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

    def set_final_sched(self, day, time, course, item1, item2):
        sched = self.schedule[day]
        sched.append([time, course, item1, item2])

    def clear_all_sched(self):
        self.schedule['M'].clear()
        self.schedule['T'].clear()
        self.schedule['W'].clear()
        self.schedule['Th'].clear()
        self.schedule['F'].clear()
        self.schedule['S'].clear()

    def clear_and_add_new_sched(self, day, sched):
        self.schedule[day] = copy.copy(sched)

    def get_schedule(self, day):
        return self.schedule[day]

    def remove_sched(self, day, value):
        self.schedule[day].remove(value)

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

        for i in range(len(population_class.get_courses())):
            class_slot = population_class.get_class_slots().pop(
                random.randrange(len(population_class.get_class_slots())))
            course = population_class.get_courses().pop(random.randrange(len(population_class.get_courses())))
            # block = population_class.get_blocks().pop(random.randrange(len(population_class.get_blocks())))
            # professor = population_class.get_professors().pop(random.randrange(len(population_class.get_professors())))

            block = course.get_block()
            professor = course.get_prof()

            gene = Gene(class_slot, course, block, professor)

            # set schedule for professor
            found_prof = data.find_prof(professor)
            if found_prof is not None:
                found_prof.set_schedule(class_slot[1], class_slot[2])
                found_prof.set_schedule_course(course.get_code(), course.get_hour())
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
                found_block.set_schedule_course(course.get_code(), course.get_hour())

            self.chromosome.add_gene(gene)

        self.class_slot_population = population_class.get_class_slots()
        GeneticAlgoPopulation.chromosome.append(self.chromosome)

    def calculate_fit(self):
        # get the current chromosome

        if len(GeneticAlgoPopulation.chromosome) == 1:
            self.curr_chromosome = GeneticAlgoPopulation.chromosome[0]
        else:
            self.curr_chromosome = GeneticAlgoPopulation.chromosome[-1]

        sum = 0
        self.chromosome_size = len(self.curr_chromosome.get_genes())

        # calculate gene fitness score
        # calculate gene fitness score
        for gene in self.curr_chromosome.get_genes():
            gene_fitness = GeneFitness(gene)
            fitness_score = self.fitness_value(gene_fitness)

            gene.add_fitness_score(fitness_score)
            sum = sum + fitness_score

        self.curr_chromosome.set_fitness_value(sum / self.chromosome_size)

        return self.curr_chromosome.get_fitness_value()

    def roulette_wheel(self):

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
            self.chromosome_size = self.chromosome_size - 1

        return selected_gene

    def hill_climbing(self, parent1):

        if parent1.get_fitness_score() == 1:
            return None, None
        else:
            # Determine which attributes can be improved
            attributes_to_improve = GeneticAlgorithm.check_attribute()
            if self.chromosome_size > 0:
                best = self.parent_two(self.curr_chromosome.get_genes(), attributes_to_improve, parent1)

                return best, attributes_to_improve
            else:
                return 0, 0

    @staticmethod
    def check_attribute():
        my_list = ['c', 'cs', 'p', 'b']

        num_elements = random.randint(1, 3)  # Randomly choose number of elements (1-3)
        search_attribute = random.sample(my_list, num_elements)  # Randomly select elements from the list

        return search_attribute

    def parent_two(self, items, attribute, parent1):

        items = copy.copy(items)

        perfect_gene_found = None
        perfect_gene_score = 0
        worst_gene_score = 1
        prev_time_slot = parent1.get_class_slot_one(2)
        prev_day = parent1.get_class_slot_one(1)

        orig_block_sched = copy.copy(parent1.get_block().get_schedule().get_schedule(prev_day))
        orig_prof_sched = copy.copy(parent1.get_professor().get_schedule().get_schedule(prev_day))
        orig_prof_course_sched = copy.copy(parent1.get_professor().get_schedule_course())
        orig_block_course_sched = copy.copy(parent1.get_block().get_schedule_course())

        size = int(self.chromosome_size * 0.6)
        count = self.chromosome_size

        for i in range(size):
            item = items.pop(random.randrange(count))

            different_day = False
            remove_prof_sched = False
            remove_block_sched = False

            parent1_copy = copy.copy(parent1)
            item_prev_day = item.get_class_slot_one(1)

            item_orig_course = copy.copy(item.get_course())
            item_orig_prof_sched1 = copy.copy(item.get_professor().get_schedule().get_schedule(item_prev_day))
            item_orig_block_sched1 = copy.copy(item.get_block().get_schedule().get_schedule(item_prev_day))
            item_orig_prof_course_sched = copy.copy(item.get_professor().get_schedule_course())
            item_orig_block_course_sched = copy.copy(item.get_block().get_schedule_course())

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

                parent1_copy.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                parent1_copy.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)

                remove_prof_sched = True
                remove_block_sched = True

                parent1_copy.get_block().set_schedule(new_day, new_time_slot)
                parent1_copy.get_professor().set_schedule(new_day, new_time_slot)

            if "c" in attribute:
                p1_course = parent1_copy.get_course()
                prev_course = parent1_copy.get_course().get_code()
                prev_course_hr = parent1_copy.get_course().get_hour()

                parent1_copy.set_course(item.get_course())
                item.set_course(p1_course)

                new_course = parent1_copy.get_course().get_code()
                new_course_hr = parent1_copy.get_course().get_hour()

                parent1_copy.get_professor().remove_schedule_course(prev_course, prev_course_hr)
                parent1_copy.get_professor().set_schedule_course(new_course, new_course_hr)
                item.get_professor().remove_schedule_course(new_course, new_course_hr)
                item.get_professor().set_schedule_course(prev_course, prev_course_hr)

                parent1_copy.get_block().remove_schedule_course(prev_course, prev_course_hr)
                parent1_copy.get_block().set_schedule_course(new_course, new_course_hr)
                item.get_block().remove_schedule_course(new_course, new_course_hr)
                item.get_block().set_schedule_course(prev_course, prev_course_hr)

            if "p" in attribute:

                item_course = item.get_course().get_code()
                item_course_hr = item.get_course().get_hour()

                parent1_copy.set_professor(item.get_professor())

                curr_course = parent1_copy.get_course().get_code()
                curr_course_hr = parent1_copy.get_course().get_hour()

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

                item.get_professor().remove_schedule_course(item_course, item_course_hr)
                parent1_copy.get_professor().set_schedule_course(curr_course, curr_course_hr)

            if "b" in attribute:

                item_course = item.get_course().get_code()
                item_course_hr = item.get_course().get_hour()

                parent1_copy.set_block(item.get_block())

                curr_course = parent1_copy.get_course().get_code()
                curr_course_hr = parent1_copy.get_course().get_hour()

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

                item.get_block().remove_schedule_course(item_course, item_course_hr)
                parent1_copy.get_block().set_schedule_course(curr_course, curr_course_hr)

            gene_fitness = GeneFitness(parent1_copy)
            fitness_score = GeneticAlgorithm.fitness_value(gene_fitness)

            parent1.get_block().get_schedule().clear_and_add_new_sched(prev_day, orig_block_sched)
            parent1.get_block().change_course_sched(orig_block_course_sched)
            parent1.get_professor().get_schedule().clear_and_add_new_sched(prev_day, orig_prof_sched)
            parent1.get_professor().change_course_sched(orig_prof_course_sched)
            item.get_block().get_schedule().clear_and_add_new_sched(item_prev_day, item_orig_block_sched1)
            item.get_block().change_course_sched(item_orig_block_course_sched)
            item.get_professor().get_schedule().clear_and_add_new_sched(item_prev_day, item_orig_prof_sched1)
            item.get_professor().change_course_sched(item_orig_prof_course_sched)
            item.set_course(item_orig_course)

            if different_day:
                parent1.get_block().get_schedule().clear_and_add_new_sched(item_prev_day, orig_block_sched2)
                parent1.get_professor().get_schedule().clear_and_add_new_sched(item_prev_day, orig_prof_sched2)
                item.get_block().get_schedule().clear_and_add_new_sched(prev_day, item_orig_block_sched2)
                item.get_professor().get_schedule().clear_and_add_new_sched(prev_day, item_orig_prof_sched2)

            if fitness_score == 1:
                return item
            elif fitness_score > perfect_gene_score:
                perfect_gene_found = item
                perfect_gene_score = fitness_score

            count -= 1

        if parent1.get_fitness_score() >= perfect_gene_score:
            return None

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
                        gene_fitness.assigned_course_to_prof() +
                        gene_fitness.assigned_course_to_block()
                ) / 12
        )

        return fitness_score

    def uniform_crossover(self, parent1, attribute, parent2):

        parent1_copy = copy.copy(parent1)
        parent2_copy = copy.copy(parent2)

        remove_prof_sched = False
        remove_block_sched = False

        if 'cs' in attribute:
            parent1.set_class_slot(parent2_copy.get_class_slot())
            parent2.set_class_slot(parent1_copy.get_class_slot())

            # get the previous time slot for parent 1 and parent 2
            p1_prev_time_slot = parent1_copy.get_class_slot_one(2)
            p1_prev_day = parent1_copy.get_class_slot_one(1)
            p2_prev_time_slot = parent2_copy.get_class_slot_one(2)
            p2_prev_day = parent2_copy.get_class_slot_one(1)

            # get the new time slot for parent 1 and parent 2
            p1_new_time_slot = parent1.get_class_slot_one(2)
            p1_new_day = parent1.get_class_slot_one(1)
            p2_new_time_slot = parent2.get_class_slot_one(2)
            p2_new_day = parent2.get_class_slot_one(1)

            # prof and block sched will change
            parent1_copy.get_block().get_schedule().remove_sched(p1_prev_day, p1_prev_time_slot)
            parent1_copy.get_professor().get_schedule().remove_sched(p1_prev_day, p1_prev_time_slot)
            parent1.get_block().set_schedule(p1_new_day, p1_new_time_slot)
            parent1.get_professor().set_schedule(p1_new_day, p1_new_time_slot)

            parent2_copy.get_block().get_schedule().remove_sched(p2_prev_day, p2_prev_time_slot)
            parent2_copy.get_professor().get_schedule().remove_sched(p2_prev_day, p2_prev_time_slot)
            parent2.get_block().set_schedule(p2_new_day, p2_new_time_slot)
            parent2.get_professor().set_schedule(p2_new_day, p2_new_time_slot)

            remove_prof_sched = True
            remove_block_sched = True

        if 'p' in attribute:

            # get the current course and alloted time for parent 1
            p1_course = parent1.get_course().get_code()
            p1_course_hour = parent1.get_course().get_hour()

            # get the current course and alloted time for parent 2
            p2_course = parent2.get_course().get_code()
            p2_course_hour = parent2.get_course().get_hour()

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

            # magpapalit ang course handle ng prof ng p1 at p2
            parent1_copy.get_professor().remove_schedule_course(p1_course, p1_course_hour)
            parent1.get_professor().set_schedule_course(p1_course, p1_course_hour)
            parent2_copy.get_professor().remove_schedule_course(p2_course, p2_course_hour)
            parent2.get_professor().set_schedule_course(p2_course, p2_course_hour)

        if 'b' in attribute:

            p1_course = parent1.get_course().get_code()
            p1_course_hour = parent1.get_course().get_hour()

            p2_course = parent2.get_course().get_code()
            p2_course_hour = parent2.get_course().get_hour()

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

            # magpapalit ang course handle ng block ng p1 at p2
            parent1_copy.get_block().remove_schedule_course(p1_course, p1_course_hour)
            parent1.get_block().set_schedule_course(p1_course, p1_course_hour)
            parent2_copy.get_block().remove_schedule_course(p2_course, p2_course_hour)
            parent2.get_block().set_schedule_course(p2_course, p2_course_hour)

        if 'c' in attribute:
            p1_course = parent1.get_course().get_code()
            p1_course_hour = parent1.get_course().get_hour()

            p2_course = parent2.get_course().get_code()
            p2_course_hour = parent2.get_course().get_hour()

            parent1.set_course(parent2_copy.get_course())
            parent2.set_course(parent1_copy.get_course())

            # magpapalit ang course handle ng block ng p1 at p2
            parent1.get_block().remove_schedule_course(p1_course, p1_course_hour)
            parent1.get_block().set_schedule_course(p2_course, p2_course_hour)
            parent2.get_block().remove_schedule_course(p2_course, p2_course_hour)
            parent2.get_block().set_schedule_course(p1_course, p1_course_hour)

            parent1.get_professor().remove_schedule_course(p1_course, p1_course_hour)
            parent1.get_professor().set_schedule_course(p2_course, p2_course_hour)
            parent2.get_professor().remove_schedule_course(p2_course, p2_course_hour)
            parent2.get_professor().set_schedule_course(p1_course, p1_course_hour)

        parent1_fitness = GeneFitness(parent1)
        parent1_fitness_score = GeneticAlgorithm.fitness_value(parent1_fitness)
        parent1.add_fitness_score(parent1_fitness_score)

        parent2_fitness = GeneFitness(parent2)
        parent2_fitness_score = GeneticAlgorithm.fitness_value(parent2_fitness)
        parent2.add_fitness_score(parent2_fitness_score)

        if parent1.get_fitness_score() == 1 and parent2.get_fitness_score() == 1:
            parent2 = self.curr_chromosome.get_genes().pop(self.curr_chromosome.get_position(parent2))
            self.add_to_new_chromosome(parent1)
            self.add_to_new_chromosome(parent2)
            self.chromosome_size = self.chromosome_size - 1
        elif parent1.get_fitness_score() == 1 and parent2.get_fitness_score() != 1:
            self.add_to_new_chromosome(parent1)
        elif parent1.get_fitness_score() != 1 and parent2.get_fitness_score() == 1:
            parent2 = self.curr_chromosome.get_genes().pop(self.curr_chromosome.get_position(parent2))
            self.add_to_new_chromosome(parent2)
            self.curr_chromosome.add_gene(parent1)
        elif parent1.get_fitness_score() > parent2.get_fitness_score():
            self.add_to_new_chromosome(parent1)
        elif parent1.get_fitness_score() < parent2.get_fitness_score():
            parent2_copy = copy.copy(parent2)
            self.curr_chromosome.pop_gene(self.curr_chromosome.get_position(parent2))
            self.curr_chromosome.add_gene(parent1)
            self.add_to_new_chromosome(parent2_copy)
        else:
            # child = self.mutation(parent1)
            self.add_to_new_chromosome(parent1)

    def prof_mutation(self):

        if self.chromosome_size != 0:

            child = random.choice(self.curr_chromosome.get_genes())
            if child.get_assigned_course_to_prof_fitness() == 0:
                course = child.get_course().get_code()
                course_hour = child.get_course().get_hour()
                prev_prof = child.get_professor()

                prev_time_slot = child.get_class_slot_one(2)
                prev_day = child.get_class_slot_one(1)

                child.set_professor(child.get_course().get_prof())

                prev_prof.get_schedule().remove_sched(prev_day, prev_time_slot)
                child.get_professor().set_schedule(prev_day, prev_time_slot)

                # change the assigned course
                prev_prof.remove_schedule_course(course, course_hour)
                child.get_professor().set_schedule_course(course, course_hour)

                fitness = GeneFitness(child)
                fitness_score = GeneticAlgorithm.fitness_value(fitness)
                child.add_fitness_score(fitness_score)

    def block_mutation(self):

        if self.chromosome_size != 0:
            child = random.choice(self.curr_chromosome.get_genes())

            if child.get_assigned_course_to_block_fitness() == 0:
                course = child.get_course().get_code()
                course_hour = child.get_course().get_hour()
                prev_block = child.get_block()

                prev_time_slot = child.get_class_slot_one(2)
                prev_day = child.get_class_slot_one(1)

                child.set_block(child.get_course().get_block())

                prev_block.get_schedule().remove_sched(prev_day, prev_time_slot)
                child.get_block().set_schedule(prev_day, prev_time_slot)

                # magpapalit ang course handle
                prev_block.remove_schedule_course(course, course_hour)
                child.get_block().set_schedule_course(course, course_hour)

                fitness = GeneFitness(child)
                fitness_score = GeneticAlgorithm.fitness_value(fitness)
                child.add_fitness_score(fitness_score)

    def mutation(self, child):

        attribute = GeneticAlgorithm.check_attribute()

        child_copy = copy.copy(child)
        prev_time_slot = child_copy.get_class_slot_one(2)
        prev_day = child_copy.get_class_slot_one(1)

        if 'cs' in attribute:

            if len(self.class_slot_population) > 0:
                random_class_slot = self.class_slot_population.pop(
                    random.randrange(len(self.class_slot_population)))

                child_class_slot = child.get_class_slot()
                self.class_slot_population.append(child_class_slot)

                child.set_class_slot(random_class_slot)

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)

                child_copy.get_class_slot_one(0).get_schedule().remove_sched(prev_day, prev_time_slot)
                child_copy.get_block().get_schedule().remove_sched(prev_day, prev_time_slot)
                child_copy.get_professor().get_schedule().remove_sched(prev_day, prev_time_slot)
                child.get_block().set_schedule(new_day, new_time_slot)
                child.get_professor().set_schedule(new_day, new_time_slot)
                child.get_class_slot_one(0).set_schedule(new_day, new_time_slot)

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
        self.prof_handled_course_sched_fitness = 0
        self.block_handled_course_sched_fitness = 0
        self.assigned_course_to_block_fitness = 0
        self.assigned_course_to_prof_fitness = 0

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

    def get_block_handled_course_sched_fitness(self):
        return self.block_handled_course_sched_fitness

    def get_prof_handled_course_sched_fitness(self):
        return self.prof_handled_course_sched_fitness

    def get_assigned_course_to_block_fitness(self):
        return self.assigned_course_to_block_fitness

    def get_assigned_course_to_prof_fitness(self):
        return self.assigned_course_to_prof_fitness

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

    def set_prof_handled_course_sched_fitness(self, value):
        self.prof_handled_course_sched_fitness = value

    def set_block_handled_course_sched_fitness(self, value):
        self.block_handled_course_sched_fitness = value

    def set_assigned_course_to_block_fitness(self, value):
        self.assigned_course_to_block_fitness = value

    def set_assigned_course_to_prof_fitness(self, value):
        self.assigned_course_to_prof_fitness = value


class GeneFitness:

    def __init__(self, gene):
        self.gene = gene

    def classroom_capacity(self):
        if self.gene.get_block().get_enrolled_students() <= self.gene.get_class_slot_one(0).get_capacity():
            self.gene.set_classroom_capacity_fitness(1)
            return 1
        else:
            self.gene.set_classroom_capacity_fitness(0)
            return 0

    def professor_work_load(self):
        day = self.gene.get_class_slot_one(1)
        schedule = self.gene.get_professor().get_schedule().get_schedule(day)
        prof_id = self.gene.get_professor().get_prof_id()

        if MEMOIZATION:
            for sublist in PWL_CACHE:
                if sublist[0] == [prof_id, day]:
                    self.gene.set_professor_work_load_fitness(sublist[1])
                    return sublist[1]

        result = has_overlap_and_multiple_copies(self.gene.get_class_slot_one(2), schedule)

        if result:
            self.gene.set_professor_work_load_fitness(0)
            if MEMOIZATION:
                PWL_CACHE.append([[prof_id, day], 0])
            return 0
        else:
            self.gene.set_professor_work_load_fitness(1)
            if MEMOIZATION:
                PWL_CACHE.append([[prof_id, day], 1])
            return 1

    def schedule_availability(self):
        day = self.gene.get_class_slot_one(1)
        room_schedule = self.gene.get_class_slot_one(0).get_schedule().get_schedule(day)
        block_schedule = self.gene.get_block().get_schedule().get_schedule(day)

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
        course = self.gene.get_course()
        available_rooms = course.get_available_rooms()

        if room.get_code() in available_rooms:
            self.gene.set_room_suitability_fitness(1)
            return 1
        else:
            self.gene.set_room_suitability_fitness(0)
            return 0

    def course_slot_suitability(self):
        slot = self.gene.get_class_slot_one(2)

        if len(slot) - 1 == self.gene.get_course().get_hour():
            self.gene.set_course_slot_suitability_fitness(1)
            return 1
        else:
            self.gene.set_course_slot_suitability_fitness(0)
            return 0

    def first_year(self):
        if self.gene.get_block().get_year() == 1:
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
        professor_schedule = self.gene.get_professor().schedule.get_schedule(prof_sched_day)
        prof_id = self.gene.get_professor().get_prof_id()

        if MEMOIZATION:
            for sublist in PL_CACHE:
                if sublist[0] == [prof_id, prof_sched_day]:
                    self.gene.set_prof_lunch_break_fitness(sublist[1])
                    return sublist[1]

        # if true then no lunch
        if has_no_lunch_break(professor_schedule):
            self.gene.set_prof_lunch_break_fitness(0)
            if MEMOIZATION:
                PL_CACHE.append([[prof_id, prof_sched_day], 0])
            return 0
        else:
            self.gene.set_prof_lunch_break_fitness(1)
            if MEMOIZATION:
                PL_CACHE.append([[prof_id, prof_sched_day], 1])
            return 1

    def block_lunch(self):
        block_code = self.gene.get_block().get_block_code()
        block_sched_day = self.gene.get_class_slot_one(1)
        block_schedule = self.gene.get_block().schedule.get_schedule(block_sched_day)

        if MEMOIZATION:
            for sublist in BL_CACHE:
                if sublist[0] == [block_code, block_sched_day]:
                    self.gene.set_block_lunch_break_fitness(sublist[1])
                    return sublist[1]

        # if true then no lunch
        if has_no_lunch_break(block_schedule):
            self.gene.set_block_lunch_break_fitness(0)
            if MEMOIZATION:
                BL_CACHE.append([[block_code, block_sched_day], 0])
            return 0
        else:
            self.gene.set_block_lunch_break_fitness(1)
            if MEMOIZATION:
                BL_CACHE.append([[block_code, block_sched_day], 1])
            return 1

    def maximum_working_hours(self):
        prof_id = self.gene.get_professor().get_prof_id()
        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        if MEMOIZATION:
            for sublist in WH_CACHE:
                if sublist[0] == [prof_id, day]:
                    self.gene.set_working_hours_fitness(sublist[1])
                    return sublist[1]

        total = 0
        for sublist in professor_schedule:
            length = len(sublist) - 1
            total += length

        if total > 8:
            total = 0
            self.gene.set_working_hours_fitness(0)
            if MEMOIZATION:
                WH_CACHE.append([[prof_id, day], 0])
        else:
            self.gene.set_working_hours_fitness(1)
            total = 1
            if MEMOIZATION:
                WH_CACHE.append([[prof_id, day], 1])

        return total

    def assigned_course_to_prof(self):
        assigned_prof = self.gene.get_course().get_prof().get_prof_id()

        if assigned_prof == self.gene.get_professor().get_prof_id():
            self.gene.set_assigned_course_to_prof_fitness(1)
            return 1
        else:
            self.gene.set_assigned_course_to_prof_fitness(0)
            return 0

    def assigned_course_to_block(self):
        assigned_block = self.gene.get_course().get_block().get_block_code()

        if assigned_block == self.gene.get_block().get_block_code():
            self.gene.set_assigned_course_to_block_fitness(1)
            return 1
        else:
            self.gene.set_assigned_course_to_block_fitness(0)
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
        self.prof = None
        self.block = None

    def get_code(self):
        return self.code

    def get_hour(self):
        return self.hour

    def get_available_rooms(self):
        return self.available_rooms

    def set_prof(self, value):
        self.prof = value

    def set_block(self, value):
        self.block = value

    def get_prof(self):
        return self.prof

    def get_block(self):
        return self.block


class LectureCourse:

    def __init__(self, code, hour, available_rooms):
        self.code = code
        self.hour = hour
        self.available_rooms = available_rooms
        self.prof = None
        self.block = None

    def get_code(self):
        return self.code

    def get_hour(self):
        return self.hour

    def get_available_rooms(self):
        return self.available_rooms

    def set_prof(self, value):
        self.prof = value

    def set_block(self, value):
        self.block = value

    def get_prof(self):
        return self.prof

    def get_block(self):
        return self.block


class Professor:

    def __init__(self, prof_id, course_code_handle):
        self.prof_id = prof_id
        self.course_code_handle = course_code_handle
        self.schedule_course = []
        self.schedule = Schedule()

    def get_prof_id(self):
        return self.prof_id

    def get_course_code_handle(self):
        new_list = [course.split("_")[-1] for course in self.course_code_handle]
        return new_list

    def get_course_block(self):
        return self.course_code_handle

    def get_schedule(self):
        return self.schedule

    def set_schedule(self, day, time_slot):
        self.schedule.set_schedule(day, time_slot)

    def set_schedule_course(self, course, hour):
        self.schedule_course.append([course, hour])

    def get_schedule_course(self):
        return self.schedule_course

    def remove_schedule_course(self, course, hour):
        self.schedule_course.remove([course, hour])

    def change_course_sched(self, sched):
        self.schedule_course = copy.copy(sched)


class Block:

    def __init__(self, block_code, enrolled_students, courses, year):
        self.block_code = block_code
        self.enrolled_students = enrolled_students
        self.courses = courses
        self.schedule_course = []
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

    def set_schedule_course(self, course, hour):
        self.schedule_course.append([course, hour])

    def get_schedule_course(self):
        return self.schedule_course

    def remove_schedule_course(self, course, hour):
        self.schedule_course.remove([course, hour])

    def change_course_sched(self, sched):
        self.schedule_course = copy.copy(sched)


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
                        c = Data.find_lec_course(course)
                        c1 = copy.copy(c)
                        coure_wo_id = block.get_block_code().replace(" ", "")
                        combine = coure_wo_id + "_" + course

                        for prof in Data.professor:
                            for prof_course in prof.get_course_block():
                                if combine == prof_course:
                                    self.courses.append(c1)
                                    self.blocks.append(block)

                                    c1.set_prof(prof)
                                    c1.set_block(block)

                for data_lab_course in Data.lab_course:
                    if course == data_lab_course.get_code():
                        courses = Data.find_lab_course(course)

                        coure_wo_id = block.get_block_code().replace(" ", "")
                        combine = coure_wo_id + "_" + course

                        for prof in Data.professor:
                            for prof_course in prof.get_course_block():
                                if combine == prof_course:
                                    c = copy.copy(courses[0])
                                    c1 = copy.copy(courses[1])

                                    self.courses.append(c)
                                    self.courses.append(c1)
                                    self.blocks.append(block)
                                    self.blocks.append(block)

                                    c.set_prof(prof)
                                    c.set_block(block)
                                    c1.set_prof(prof)
                                    c1.set_block(block)

                        break

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


def convert_time(start_time, end_time):
    start_time = int(start_time)
    end_time = int(end_time)

    start_hour = start_time % 12
    end_hour = end_time % 12

    if start_hour == 0:
        start_hour = 12

    if end_hour == 0:
        end_hour = 12

    start_suffix = 'pm' if start_time >= 12 else 'am'
    end_suffix = 'pm' if end_time >= 12 else 'am'

    start_str = '{}:00 {}'.format(start_hour, start_suffix)
    end_str = '{}:00 {}'.format(end_hour, end_suffix)

    return '{} to {}'.format(start_str, end_str)


def reformat_sched(lst):
    new_lst = []
    separator = "*"
    count = 1

    if len(lst) == 1:
        separator = ""

    for sched in lst:

        if count == len(lst):
            separator = ""

        new_time = convert_time(sched[0][0], sched[0][-1])
        format = new_time + "\n" + sched[1] + " - " + sched[2] + " (" + sched[3] + ")" + "\n " + separator
        new_lst.append(format)

        count += 1

    return new_lst


def display_table(lst, obj):
    for i in range(len(lst)):

        sched_row = []
        data = []

        if obj == "PROFESSOR":
            id = lst[i].get_prof_id()
        elif obj == "BLOCK":
            id = lst[i].get_block_code()
        else:
            id = lst[i].get_code()

        monday_sched = sorted(lst[i].get_schedule().get_m(), key=lambda x: int(x[0][0]))
        monday_sched = reformat_sched(monday_sched)

        tuesday_sched = sorted(lst[i].get_schedule().get_t(), key=lambda x: int(x[0][0]))
        tuesday_sched = reformat_sched(tuesday_sched)

        wednesday_sched = sorted(lst[i].get_schedule().get_w(), key=lambda x: int(x[0][0]))
        wednesday_sched = reformat_sched(wednesday_sched)

        thursday_sched = sorted(lst[i].get_schedule().get_th(), key=lambda x: int(x[0][0]))
        thursday_sched = reformat_sched(thursday_sched)

        friday_sched = sorted(lst[i].get_schedule().get_f(), key=lambda x: int(x[0][0]))
        friday_sched = reformat_sched(friday_sched)

        saturday_sched = sorted(lst[i].get_schedule().get_s(), key=lambda x: int(x[0][0]))
        saturday_sched = reformat_sched(saturday_sched)

        list_lengths = [len(monday_sched), len(tuesday_sched), len(wednesday_sched), len(thursday_sched),
                        len(friday_sched), len(saturday_sched)]
        greatest_lenth = max(list_lengths)

        print(id + " SCHEDULE")
        for row in range(greatest_lenth):
            row_list = []
            row_list.append("" if row >= len(monday_sched) else monday_sched[row])
            row_list.append("" if row >= len(tuesday_sched) else tuesday_sched[row])
            row_list.append("" if row >= len(wednesday_sched) else wednesday_sched[row])
            row_list.append("" if row >= len(thursday_sched) else thursday_sched[row])
            row_list.append("" if row >= len(friday_sched) else friday_sched[row])
            row_list.append("" if row >= len(saturday_sched) else saturday_sched[row])
            sched_row.append(row_list)

        for row in sched_row:
            data.append(row)

        print(tabulate(data, headers=["MONDAY", "TUESDAY", "WEDNESDAY", "THURSDAY", "FRIDAY", "SATURDAY"],
                       tablefmt="pretty"))
        print("\n\n")


if __name__ == '__main__':
    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F"]

    sl1 = ClassRoom("SL1", 50, "Laboratory")
    sl2 = ClassRoom("SL2", 50, "Laboratory")
    sl3 = ClassRoom("SL3", 50, "Laboratory")
    cisco = ClassRoom("CISCO", 50, "Laboratory")
    ml = ClassRoom("ML", 50, "Laboratory")
    itl = ClassRoom("ITL", 50, "Laboratory")
    edl = ClassRoom("EDL", 50, "Laboratory")
    r102 = ClassRoom("102", 50, "Regular")
    r103 = ClassRoom("103", 50, "Regular")
    r104 = ClassRoom("104", 50, "Regular")
    r201 = ClassRoom("201", 50, "Regular")
    r202 = ClassRoom("202", 50, "Regular")
    r203 = ClassRoom("203", 50, "Regular")
    r401 = ClassRoom("401", 50, "Regular")
    r402 = ClassRoom("402", 50, "Regular")
    r403 = ClassRoom("403", 50, "Regular")

    Data.add_room(sl1)
    Data.add_room(sl2)
    Data.add_room(sl3)
    Data.add_room(cisco)
    Data.add_room(ml)
    Data.add_room(itl)
    Data.add_room(edl)
    Data.add_room(r102)
    Data.add_room(r103)
    Data.add_room(r104)
    Data.add_room(r201)
    Data.add_room(r202)
    Data.add_room(r203)
    Data.add_room(r401)
    Data.add_room(r402)
    Data.add_room(r403)

    Data.add_lab_course(LabCourse("CS 121", 3, ["SL1", "SL2", "SL3", "ML"]))
    Data.add_lab_course(LabCourse("CS 111", 3, ["SL1", "SL2", "SL3", "ML"]))
    Data.add_lab_course(LabCourse("CS 131", 3, ["SL1", "SL2", "SL3", "ML", "ITL"]))

    Data.add_lab_course(
        LabCourse("CS 121", 2, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lab_course(
        LabCourse("CS 111", 2, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lab_course(
        LabCourse("CS 131", 2, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))

    Data.add_lec_course(
        LectureCourse("FILI 102", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lec_course(
        LectureCourse("GEd 105", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lec_course(
        LectureCourse("GEd 108", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lec_course(
        LectureCourse("MATH 401", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lec_course(
        LectureCourse("MATH 407", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lec_course(
        LectureCourse("MATH 111", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))
    Data.add_lec_course(
        LectureCourse("GEd 109", 3, ["101", "102", "103", "104", "201", "202", "203", "EDL", "401", "402", "403"]))

    cs1201 = Block("CS 1201", 40, ["CS 121", "FILI 102", "GEd 105", "GEd 108", "MATH 401", "MATH 407"], 1)
    cs1202 = Block("CS 1202", 45, ["CS 121", "FILI 102", "GEd 105", "GEd 108", "MATH 401", "MATH 407"], 1)
    cs1203 = Block("CS 1203", 43, ["CS 121", "FILI 102", "GEd 105", "GEd 108", "MATH 401", "MATH 407"], 1)
    it1201 = Block("IT 1201", 50, ["CS 111", "CS 131", "MATH 111", "FILI 102", "GEd 105", "GEd 109"], 1)
    it1202 = Block("IT 1202", 45, ["CS 111", "CS 131", "MATH 111", "FILI 102", "GEd 105", "GEd 109"], 1)
    it1203 = Block("IT 1203", 40, ["CS 111", "CS 131", "MATH 111", "FILI 102", "GEd 105", "GEd 109"], 1)
    it1204 = Block("IT 1204", 45, ["CS 111", "CS 131", "MATH 111", "FILI 102", "GEd 105", "GEd 109"], 1)
    it1205 = Block("IT 1205", 45, ["CS 111", "CS 131", "MATH 111", "FILI 102", "GEd 105", "GEd 109"], 1)
    it1206 = Block("IT 1206", 40, ["CS 111", "CS 131", "MATH 111", "FILI 102", "GEd 105", "GEd 109"], 1)

    Data.add_block(cs1201)
    Data.add_block(cs1202)
    Data.add_block(cs1203)
    Data.add_block(it1201)
    Data.add_block(it1202)
    Data.add_block(it1203)
    Data.add_block(it1204)
    Data.add_block(it1205)
    Data.add_block(it1206)

    p1 = Professor("P1", ["CS1201_CS 121", "CS1202_CS 121", "CS1201_MATH 401", "IT1201_CS 111", "IT1201_GEd 109"])
    p2 = Professor("P2", ["CS1203_CS 121", "CS1201_FILI 102", "CS1201_GEd 108", "IT1206_MATH 111", "IT1201_MATH 111"])
    p3 = Professor("P3", ["CS1202_FILI 102", "CS1203_MATH 401", "IT1202_CS 111", "IT1203_CS 111", "IT1202_GEd 109"])
    p4 = Professor("P4", ["IT1201_FILI 102", "CS1201_MATH 407", "CS1202_MATH 407", "IT1202_MATH 111", "IT1205_GEd 109",
                          "IT1206_GEd 109"])
    p5 = Professor("P5", ["CS1201_GEd 105", "IT1201_CS 131", "IT1202_CS 131", "IT1203_FILI 102", "CS1202_GEd 105",
                          "CS1203_GEd 105"])
    p6 = Professor("P6", ["IT1201_GEd 105", "IT1202_GEd 105", "IT1203_CS 131", "IT1205_CS 131", "IT1203_MATH 111",
                          "IT1203_GEd 105"])
    p7 = Professor("P7", ["IT1205_CS 111", "IT1206_CS 131", "CS1203_FILI 102", "IT1202_FILI 102", "IT1203_GEd 109"])
    p8 = Professor("P8", ["CS1202_GEd 108", "CS1203_GEd 108", "IT1204_CS 111", "IT1204_FILI 102", "IT1204_GEd 109"])
    p9 = Professor("P9", ["CS1202_MATH 401", "IT1204_CS 131", "IT1204_MATH 111", "IT1205_FILI 102", "IT1206_FILI 102"])
    p10 = Professor("P10", ["CS1203_MATH 407", "IT1206_CS 111", "IT1205_MATH 111", "IT1204_GEd 105", "IT1205_GEd 105",
                            "IT1206_GEd 105"])

    Data.add_prof(p1)
    Data.add_prof(p2)
    Data.add_prof(p3)
    Data.add_prof(p4)
    Data.add_prof(p5)
    Data.add_prof(p6)
    Data.add_prof(p7)
    Data.add_prof(p8)
    Data.add_prof(p9)
    Data.add_prof(p10)

    # Some code here
    start_time = time.time()
    population = Population()
    population.initialize_population(meetingTime, day)
    end_time = time.time()

    elapsed_time = end_time - start_time
    print("---------------------------------------------------------------------")
    print("Class-slot population size : ", len(population.get_class_slots()))
    print("Professor population size : ", len(population.get_professors()))
    print("Block population size : ", len(population.get_blocks()))
    print("Course population size : ", len(population.get_courses()))
    print("---------------------------------------------------------------------")
    print(f"Population initialization elapsed time: {elapsed_time} seconds")

    start_time = time.time()
    genetic_algo = GeneticAlgorithm()
    genetic_algo.encode_chromosome(population)
    end_time = time.time()

    print(f"Encoding Chromosome elapsed time: {elapsed_time} seconds")

    count = 0
    not_perfect_schedule = True
    generation = 0

    start_time = time.time()
    while not_perfect_schedule:

        chromosome_fit = genetic_algo.calculate_fit()
        new_chromosome = True
        MEMOIZATION = False
        generation = generation + 1

        # print("Generation " + str(count) + " fitness score : ", genetic_algo.curr_chromosome.get_fitness_value())

        count = count + 1
        if chromosome_fit == 1:
            end_time = time.time()
            print(f"Finding Best Schedule elapsed time: {elapsed_time} seconds")

            print("---------------------------------------------------------------------")
            print("Solution found at generation", generation-1)
            print("Fitness score : ", genetic_algo.curr_chromosome.get_fitness_value())
            print("---------------------------------------------------------------------\n\n")

            prof_list = [p1, p2, p3, p4, p5, p6, p7, p8, p9, p10]
            room_list = [sl1, sl2, sl3, cisco, ml, itl, edl, r102, r103, r104, r201, r202, r203, r401, r402, r403]
            block_list = [cs1201, cs1202, cs1203, it1201, it1202, it1203, it1204, it1205, it1206]

            if len(GeneticAlgoPopulation.chromosome) == 1:
                genetic_algo.curr_chromosome = GeneticAlgoPopulation.chromosome[0]
            else:
                genetic_algo.curr_chromosome = GeneticAlgoPopulation.chromosome[-1]

            # delete all sched

            for professor in prof_list:
                professor.get_schedule().clear_all_sched()
            for r in room_list:
                r.get_schedule().clear_all_sched()
            for b in block_list:
                b.get_schedule().clear_all_sched()

            # set new sched from genes found
            for gene in genetic_algo.curr_chromosome.get_genes():
                day = gene.get_class_slot_one(1)
                time_slot = gene.get_class_slot_one(2)
                room = gene.get_class_slot_one(0).get_code()
                prof = gene.get_professor().get_prof_id()
                block = gene.get_block().get_block_code()
                course = gene.get_course().get_code()

                gene.get_block().get_schedule().set_final_sched(day, time_slot, course, prof, room)
                gene.get_professor().get_schedule().set_final_sched(day, time_slot, course, block, room)
                gene.get_class_slot_one(0).get_schedule().set_final_sched(day, time_slot, course, block, prof)

            print("PROFESSORS SCHEDULE\n")
            display_table(prof_list, "PROFESSOR")
            print("BLOCKS SCHEDULE\n")
            display_table(block_list, "BLOCK")
            print("ROOM SCHEDULE\n")
            display_table(room_list, "ROOM")

            not_perfect_schedule = False

        else:
            while new_chromosome:
                parent1 = genetic_algo.roulette_wheel()

                if parent1 is not None:

                    parent2, attribute = genetic_algo.hill_climbing(parent1)
                    if parent2 is None or attribute is None:
                        genetic_algo.add_to_new_chromosome(parent1)
                    elif parent2 == 0 and attribute == 0:
                        child = genetic_algo.mutation(parent1)
                        genetic_algo.add_to_new_chromosome(child)
                    else:
                        genetic_algo.uniform_crossover(parent1, attribute, parent2)

                        genetic_algo.block_mutation()
                        genetic_algo.prof_mutation()
                else:
                    print("PARENT 1 IS NONE")

                if genetic_algo.get_chromosome_size() == 0:
                    genetic_algo.add_new_chromosome()
                    WH_CACHE.clear()
                    PL_CACHE.clear()
                    BL_CACHE.clear()
                    PWL_CACHE.clear()

                    MEMOIZATION = True
                    new_chromosome = False
