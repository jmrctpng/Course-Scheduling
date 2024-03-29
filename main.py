from tabulate import tabulate
from typing import List
import random
import time
import numpy as np
import copy
import pandas as pd

PROF_MUTATE = None
PROF_MUTATE1 = None
BLOCK_MUTATE = None
BLOCK_MUTATE1 = None
WH_CACHE = []
PL_CACHE = []
BL_CACHE = []
PWL_CACHE = []
MEMOIZATION = True

PROF_LIST = []
BLOCK_LIST = []
ROOM_LIST = []

LAB_COURSE = []
REG_COURSE = []


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

            choices = [0,1]
            probabilities = [0.60, 0.40]

            res = random.choices(choices, probabilities)

            if res[0] == 0:
                if self.chromosome_size > 0:
                    best = self.parent_two(self.curr_chromosome.get_genes(), attributes_to_improve, parent1)
                    return best, attributes_to_improve
                else:
                    return 0, 0
            else:
                return 0,0

    @staticmethod
    def check_attribute():
        return ['cs']

    def parent_two(self, items, attribute, parent1):

        items = copy.copy(items)

        perfect_gene_found = None
        perfect_gene_score = 0
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

                parent1_copy.get_block().set_schedule(new_day, new_time_slot)
                parent1_copy.get_professor().set_schedule(new_day, new_time_slot)

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
            self.add_to_new_chromosome(parent1)

    def mutation(self, child):

        found = False
        curr_fitness = child.get_fitness_score()
        size = len(self.class_slot_population)
        count = 0

        if len(self.class_slot_population) > 0:
            while not found:

                child_copy = copy.copy(child)
                prev_time_slot = child_copy.get_class_slot_one(2)
                prev_day = child_copy.get_class_slot_one(1)

                random_class_slot = self.class_slot_population.pop(
                    random.randrange(len(self.class_slot_population)))

                child_class_slot = child_copy.get_class_slot()
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

                gene_fitness = GeneFitness(child)
                fitness_score = GeneticAlgorithm.fitness_value(gene_fitness)

                if fitness_score == 1 or fitness_score > curr_fitness:
                    found = True

                if count == size:
                    found = True

                count += 1

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

    def add_all_genes(self, genes):
        self.genes = genes

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
        self.copy_courses = copy.copy(courses)
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

    def remove_course(self, course):
        self.copy_courses.remove(course)

    def get_copy_courses(self):
        return self.copy_courses


class Population:
    def __init__(self):
        self.population = []
        self.courses = []
        self.blocks = []
        self.professors = []
        self.class_slots = []

    def initialize_population(self, meetingTime, meetingDay):
        self.generate_prof_population()
        self.generate_course_block_population()
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


def assign_section(course_code):
    blocks = Data.block

    for block in blocks:
        for course in block.get_copy_courses():
            if course == course_code:
                block.remove_course(course)
                id = block.get_block_code().upper().replace(" ", "") + "_" + course
                return id


def generate_lab_subjects(df):
    lab_subj = []

    # iterate through each row and column to access cell values
    for index, row in df.iterrows():

        i = 0
        id = ""
        hours = 0
        rooms = []

        for col in df.columns:
            cell = row[col]
            if not pd.isna(cell):
                if i == 0:
                    id = cell
                elif i == 1:
                    hours = cell
                else:
                    rooms.append(cell)
                i += 1

        lab_subj.append([id, hours, rooms])

    for subj in lab_subj:
        subj = LabCourse(subj[0], subj[1], subj[2])
        Data.add_lab_course(subj)
        LAB_COURSE.append(subj.get_code())


def generate_class_room(df):
    rooms = []

    # iterate through each row and column to access cell values
    for index, row in df.iterrows():

        i = 0
        id = ""
        capacity = 0
        type = ""

        for col in df.columns:
            cell = row[col]
            if not pd.isna(cell):
                if i == 0:
                    id = cell
                elif i == 1:
                    capacity = cell
                elif i == 2:
                    type = cell
                i += 1

        rooms.append([id, capacity, type])

    for room in rooms:
        r = ClassRoom(room[0], room[1], room[2])
        Data.add_room(r)
        ROOM_LIST.append(r)


def generate_regular_subjects(df):
    reg_subj = []

    # iterate through each row and column to access cell values
    for index, row in df.iterrows():

        i = 0
        id = ""
        hours = 0
        rooms = []

        for col in df.columns:
            cell = row[col]
            if not pd.isna(cell):
                if i == 0:
                    id = cell
                elif i == 1:
                    hours = cell
                else:
                    rooms.append(cell)
                i += 1

        reg_subj.append([id, hours, rooms])

    for subj in reg_subj:
        subj = LectureCourse(subj[0], subj[1], subj[2])
        Data.add_lec_course(subj)
        REG_COURSE.append(subj.get_code())


def generate_block_list(df):
    block_list = []

    # iterate through each row and column to access cell values
    for index, row in df.iterrows():

        i = 0
        id = ""
        enrolled_students = 0
        year = 0
        subjects = []

        for col in df.columns:
            cell = row[col]
            if not pd.isna(cell):
                if i == 0:
                    id = cell
                elif i == 1:
                    enrolled_students = cell
                elif i == 2:
                    year = cell
                else:
                    subjects.append(cell)
                i += 1
        block_list.append([id, enrolled_students, subjects, year])

    for block in block_list:
        block = Block(block[0], block[1], block[2], block[3])
        Data.add_block(block)
        BLOCK_LIST.append(block)


def generate_professor_list(df):
    prof_list = []

    # iterate through each row and column to access cell values
    for index, row in df.iterrows():

        i = 0
        id = ""
        subjects = []

        for col in df.columns:
            cell = row[col]
            if not pd.isna(cell):
                if i == 0:
                    id = cell
                else:
                    course_with_assigned_block = assign_section(cell)
                    subjects.append(course_with_assigned_block)
                i += 1
        prof_list.append([id, subjects])

    for prof in prof_list:
        prof = Professor(prof[0], prof[1])
        Data.add_prof(prof)
        PROF_LIST.append(prof)


if __name__ == '__main__':
    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F", "S"]

    # read the csv file
    room_df = pd.read_csv('rooms.csv')
    generate_class_room(room_df)

    # read the csv file
    lab_df = pd.read_csv('Dataset 3/D3lab_course.csv')
    generate_lab_subjects(lab_df)

    # read the csv file
    reg_df = pd.read_csv('Dataset 3/D3reg_course.csv')
    generate_regular_subjects(reg_df)

    # read the csv file
    block_df = pd.read_csv('Dataset 3/D3block_data.csv')
    generate_block_list(block_df)

    # read the csv file
    prof_df = pd.read_csv('Dataset 3/dataset3_professors.csv')
    generate_professor_list(prof_df)

    population = Population()
    population.initialize_population(meetingTime, day)
    genetic_algo = GeneticAlgorithm()
    genetic_algo.encode_chromosome(population)

    count = 0
    best_generation_fitness_score = 0
    best_generation_index = 0
    best = None
    not_perfect_schedule = True
    generation = 0

    while not_perfect_schedule:
        chromosome_fit = genetic_algo.calculate_fit()

        new_chromosome = True
        MEMOIZATION = False
        generation = generation + 1

        print(str(chromosome_fit) + ",")

        if chromosome_fit > best_generation_fitness_score:
            best = copy.copy(genetic_algo.curr_chromosome.get_genes())
            best_generation_fitness_score = chromosome_fit
            best_generation_index = generation

        count = count + 1
        if chromosome_fit == 1:

            print("fitness score")

            best_chromosome = Chromosome()
            best_chromosome.add_all_genes(best)

            print("---------------------------------------------------------------------")
            print("Solution found at generation", best_generation_index)
            print("Fitness score : ", best_generation_fitness_score)
            print("---------------------------------------------------------------------\n\n")

            # delete all sched
            for professor in PROF_LIST:
                professor.get_schedule().clear_all_sched()
            for r in ROOM_LIST:
                r.get_schedule().clear_all_sched()
            for b in BLOCK_LIST:
                b.get_schedule().clear_all_sched()

            # set new sched from genes found
            for gene in best_chromosome.get_genes():
                day = gene.get_class_slot_one(1)
                time_slot = gene.get_class_slot_one(2)
                room = gene.get_class_slot_one(0).get_code()
                prof = gene.get_professor().get_prof_id()
                block = gene.get_block().get_block_code()
                course = gene.get_course().get_code()

                if gene.get_first_year_fitness() == 2:
                    first_year = 1
                else:
                    first_year = gene.get_first_year_fitness()

                print("Gene ID : ", gene)
                print("Class slot: ", gene.get_class_slot_one(0).get_code(), " ",
                      gene.get_class_slot_one(0).get_capacity(), " ",
                      gene.get_class_slot_one(0).get_type_of_room(), " ", gene.get_class_slot_one(1), " ",
                      gene.get_class_slot_one(2))
                print("Course code: ", gene.get_course().get_code())
                print("assigned prof: ", gene.get_course().get_prof().get_prof_id())
                print("assigned block: ", gene.get_course().get_block().get_block_code())
                print("Course hr : ", gene.get_course().get_hour())
                print("Professor: ", gene.get_professor().get_prof_id())
                print("Block: ", gene.get_block().get_block_code())
                print("Score :", gene.get_fitness_score())
                print("CR :", gene.get_classroom_capacity_fitness())
                print("PW :", gene.get_professor_work_load_fitness())
                print("RA :", gene.get_room_availability_fitness())
                print("FF :", first_year)
                print("BA :", gene.get_block_availability_fitness())
                print("RS :", gene.get_room_suitability_fitness())
                print("CSS :", gene.get_course_slot_suitability_fitness())
                print("PL :", gene.get_prof_lunch_break_fitness())
                print("BL :", gene.get_block_lunch_break_fitness())
                print("WH :", gene.get_working_hours_fitness())
                print("ACP :", gene.get_assigned_course_to_prof_fitness())
                print("ACB :", gene.get_assigned_course_to_block_fitness())
                print("=======================================")

                gene.get_block().get_schedule().set_final_sched(day, time_slot, course, prof, room)
                gene.get_professor().get_schedule().set_final_sched(day, time_slot, course, block, room)
                gene.get_class_slot_one(0).get_schedule().set_final_sched(day, time_slot, course, block, prof)

            print("PROFESSORS SCHEDULE\n")
            display_table(PROF_LIST, "PROFESSOR")
            print("BLOCKS SCHEDULE\n")
            display_table(BLOCK_LIST, "BLOCK")
            print("ROOM SCHEDULE\n")
            display_table(ROOM_LIST, "ROOM")

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
