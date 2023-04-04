from typing import List
import random
import time
import numpy as np
import copy
import threading
from collections import Counter

WH_CACHE = []
PL_CACHE = []
BL_CACHE = []
PWL_CACHE = []
BHCS_CACHE = []
PHCS_CACHE = []
MEMOIZATION = True


LAB_COURSE = ['CS 121', 'CS 222', 'CS 322', 'cs 221']

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
            print("PHCS :", gene.get_prof_handled_course_sched_fitness())
            print("BHCS :", gene.get_block_handled_course_sched_fitness())
            print("=======================================")
            print("prof id : ", gene.professor)
            print("prof sched : ", gene.professor.get_schedule().get_schedule(gene.class_slot[1]))
            print("block sched : ", gene.block.get_schedule().get_schedule(gene.class_slot[1]))
            print("Room sched : ", gene.class_slot[0].get_schedule().get_schedule(gene.class_slot[1]))
            print("-----------------------------")


        self.curr_chromosome.set_fitness_value(sum / self.chromosome_size)

        print("chromosome total fit : ", self.curr_chromosome.get_fitness_value())

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
            # self.curr_chromosome.pop_gene(index)
            self.chromosome_size = self.chromosome_size - 1
            # print("update size : ", self.chromosome_size)

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

        check_attribute = ['cs', 'b', 'c', 'p']
        search_attribute = []
        class_slot = [
            ['PW', 0 if parent1.get_professor_work_load_fitness() == 0 else 1],
            ['RA', 0 if parent1.get_room_availability_fitness() == 0 else 1],
            ['BA', 0 if parent1.get_block_availability_fitness() == 0 else 1],
            ['CSS', 0 if parent1.get_course_slot_suitability_fitness() == 0 else 1],
            ['PL', 0 if parent1.get_prof_lunch_break_fitness() == 0 else 1],
            ['BL', 0 if parent1.get_block_lunch_break_fitness() == 0 else 1],
            ['WH', 0 if parent1.get_working_hours_fitness() == 0 else 1]
        ]

        all_fitness_score = [
            ['CR', 0 if parent1.get_classroom_capacity_fitness() == 0 else 1],
            ['PW', 0 if parent1.get_professor_work_load_fitness() == 0 else 1],
            ['RA', 0 if parent1.get_room_availability_fitness() == 0 else 1],
            ['BA', 0 if parent1.get_block_availability_fitness() == 0 else 1],
            ['RS', 0 if parent1.get_room_suitability_fitness() == 0 else 1],
            ['CSS', 0 if parent1.get_course_slot_suitability_fitness() == 0 else 1],
            ['PL', 0 if parent1.get_prof_lunch_break_fitness() == 0 else 1],
            ['BL', 0 if parent1.get_block_lunch_break_fitness() == 0 else 1],
            ['WH', 0 if parent1.get_working_hours_fitness() == 0 else 1],
            ['PHC', 0 if parent1.get_prof_handled_course_fitness() == 0 else 1],
            ['BEC', 0 if parent1.get_block_enrolled_course_fitness() == 0 else 1],
            ['BHCS', 0 if parent1.get_block_handled_course_sched_fitness() == 0 else 1],
            ['PHCS', 0 if parent1.get_prof_handled_course_sched_fitness() == 0 else 1]
        ]

        courses_fitness_score = [
            ['CSS', 0 if parent1.get_course_slot_suitability_fitness() == 0 else 1],
            ['PHC', 0 if parent1.get_prof_handled_course_fitness() == 0 else 1],
            ['RS', 0 if parent1.get_room_suitability_fitness() == 0 else 1],
            ['BEC', 0 if parent1.get_block_enrolled_course_fitness() == 0 else 1],
            ['BHCS', 0 if parent1.get_block_handled_course_sched_fitness() == 0 else 1],
            ['PHCS', 0 if parent1.get_prof_handled_course_sched_fitness() == 0 else 1]
        ]

        prof_fitness_score = [
            ['PW', 0 if parent1.get_professor_work_load_fitness() == 0 else 1],
            ['PL', 0 if parent1.get_prof_lunch_break_fitness() == 0 else 1],
            ['PHC', 0 if parent1.get_prof_handled_course_fitness() == 0 else 1],
            ['PHCS', 0 if parent1.get_prof_handled_course_sched_fitness() == 0 else 1]
        ]

        block_fitness_score = [
            ['BA', 0 if parent1.get_block_availability_fitness() == 0 else 1],
            ['BL', 0 if parent1.get_block_lunch_break_fitness() == 0 else 1],
            ['BEC', 0 if parent1.get_block_enrolled_course_fitness() == 0 else 1],
            ['BHCS', 0 if parent1.get_block_handled_course_sched_fitness() == 0 else 1]
        ]

        zeros = 0
        ones = 0
        all_zeros = 0
        all_ones = 0
        courses_zero = 0
        courses_one = 0
        prof_zeros = 0
        prof_ones = 0
        block_zeros = 0
        block_ones = 0

        # class slot
        for sublist in class_slot:
            if sublist[1] == 0:
                zeros += 1
            elif sublist[1] == 1:
                ones += 1

        # all fitness
        for sublist in all_fitness_score:
            if sublist[1] == 0:
                all_zeros += 1
            elif sublist[1] == 1:
                all_ones += 1

        # courses
        for sublist in courses_fitness_score:
            if sublist[1] == 0:
                courses_zero += 1
            elif sublist[1] == 1:
                courses_one += 1

        # prof
        for sublist in prof_fitness_score:
            if sublist[1] == 0:
                prof_zeros += 1
            elif sublist[1] == 1:
                prof_ones += 1

        # block
        for sublist in block_fitness_score:
            if sublist[1] == 0:
                block_zeros += 1
            elif sublist[1] == 1:
                block_ones += 1

        block_fitness_with_zero = [element for element in block_fitness_score if element[1] == 0]

        if block_zeros >= block_ones:
            search_attribute.append("b")  # change block
        else:
            attribute_to_change = [['BA', 'cs'], ['BL', 'cs'], ['BEC', 'c'], ['BHCS', 'c']]

            for element in block_fitness_with_zero:
                for attribute in attribute_to_change:
                    if element[0] == attribute[0]:
                        search_attribute.append(attribute[1])

        courses_fitness_with_zero = [element for element in courses_fitness_score if element[1] == 0]
        if courses_zero >= courses_one:
            search_attribute.append("c")  # change course
        else:
            courses_with_zero = [element for element in courses_fitness_score if element[1] == 0]
            attribute_to_change = [['PHC', 'p'], ['BEC', 'b'], ['RS', 'r'], ['BHCS', 'b'], ['CSS', 'cs'], ['PHCS', 'p']]

            for element in courses_with_zero:
                for attribute in attribute_to_change:
                    if element[0] == attribute[0]:
                        search_attribute.append(attribute[1])

        prof_fitness_with_zero = [element for element in prof_fitness_score if element[1] == 0]

        if prof_zeros >= prof_ones:
            search_attribute.append("p")  # change course
        else:
            attribute_to_change = [['PW', 'cs'], ['PL', 'cs'], ['PHC', 'c'], ['PHCS', 'c']]

            for element in prof_fitness_with_zero:
                for attribute in attribute_to_change:
                    if element[0] == attribute[0]:
                        search_attribute.append(attribute[1])

        all_fitness_with_zero = [element for element in all_fitness_score if element[1] == 0]

        if all_ones > all_zeros:
            if len(all_fitness_with_zero) == 2:
                if all_fitness_with_zero[0][0] == 'PHC' and all_fitness_with_zero[1][0] == 'BHCS':
                    search_attribute.append('b')
                if all_fitness_with_zero[0][0] == 'BHCS' and all_fitness_with_zero[1][0] == 'PHC':
                    search_attribute.append('b')
            if len(all_fitness_with_zero) == 1:
                if all_fitness_with_zero[0][0] == 'BHCS':
                    search_attribute.append('b')

        cs_fitness_with_zero = [element for element in class_slot if element[1] == 0]

        if zeros > ones:
            search_attribute.append("cs")  # change course slot
        else:
            cs_to_change = [['PW', 'p'], ['RA', 'r'], ['BA', 'b'], ['CSS', 'c'], ['PL', 'p'], ['BL', 'b'], ['WH', 'd']]

            for element in cs_fitness_with_zero:
                for attribute in cs_to_change:
                    if element[0] == attribute[0]:
                        search_attribute.append(attribute[1])

        if ones >= zeros and courses_one >= courses_zero:
            if 'CSS' in cs_fitness_with_zero and 'CSS' in courses_fitness_with_zero:
                search_attribute.remove('c')
        if courses_one >= courses_zero and prof_ones >= prof_zeros:
            if 'PHC' in courses_fitness_with_zero and 'PHC' in prof_fitness_with_zero:
                search_attribute.remove('c')
            if 'PHCS' in courses_fitness_with_zero and 'PHCS' in prof_fitness_with_zero:
                search_attribute.remove('c')
        if ones >= zeros and prof_ones >= prof_zeros:
            if 'PW' in cs_fitness_with_zero and 'PW' in cs_fitness_with_zero:
                search_attribute.remove('cs')
            if 'PL' in cs_fitness_with_zero and 'PL' in cs_fitness_with_zero:
                search_attribute.remove('cs')
        if courses_one >= courses_zero and block_ones >= block_zeros:
            if 'BEC' in courses_fitness_with_zero and 'BEC' in block_fitness_with_zero:
                search_attribute.remove('c')
            if 'BHCS' in courses_fitness_with_zero and 'BEC' in block_fitness_with_zero:
                search_attribute.remove('c')
        if ones >= zeros and block_ones >= block_zeros:
            if 'BA' in cs_fitness_with_zero and 'BA' in block_fitness_with_zero:
                search_attribute.remove('cs')
            if 'BL' in cs_fitness_with_zero and 'BL' in cs_fitness_with_zero:
                search_attribute.remove('cs')

        if (
                parent1.get_prof_handled_course_fitness() == 1 and parent1.get_block_enrolled_course_fitness() == 1 and parent1.get_room_suitability_fitness() == 0) or (
                parent1.get_classroom_capacity_fitness() == 0):
            search_attribute.append("r")  # change room
        if parent1.get_first_year_fitness() == 0 or parent1.get_working_hours_fitness() == 0:
            search_attribute.append("d")  # change day
        if parent1.get_block_handled_course_sched_fitness() == 0 and parent1.get_block_enrolled_course_fitness() == 1:

            three_hours = 0
            two_hours = 0
            laboratory = False

            course = parent1.get_course().get_code()

            if course in LAB_COURSE:
                three_hours = 3
                two_hours = 2
                laboratory = True
            else:
                three_hours = 3
                two_hours = 0

            # get all the schedules from monday - saturday
            schedule_course_code = parent1.get_block().get_schedule_course()

            for sched in schedule_course_code:
                if sched[0] == course:
                    if sched[1] == 3:
                        three_hours -= 3
                    elif sched[1] == 2:
                        two_hours -= 2

                    if three_hours < 0 and two_hours < 0:
                        break

            if laboratory:
                if three_hours == 0 and two_hours > 0:
                    # look for same course with 2 hours time slot
                    search_attribute.append("cs")  # change course slot
                elif three_hours == 0 and two_hours < 0 and parent1.get_course().get_hour() == 2:
                    search_attribute.append("c")  # change course
                elif three_hours > 0 and two_hours == 0:
                    # look for same course with 3 hours time slot
                    search_attribute.append("cs")  # change course slot
                elif three_hours < 0 and two_hours == 0 and parent1.get_course().get_hour() == 3:
                    search_attribute.append("c")  # change course
                elif three_hours < 0 and two_hours < 0:
                    search_attribute.append("c")  # change course
                elif three_hours == 0 and two_hours == 0:
                    if parent1.get_prof_handled_course_fitness() == 0:
                        search_attribute.append("p")  # change prof
                    if parent1.get_room_suitability_fitness() == 0:
                        search_attribute.append("r")  # change room
            else:
                if three_hours == 0:
                    if parent1.get_prof_handled_course_fitness() == 0:
                        search_attribute.append("p")  # change prof
                    if parent1.get_room_suitability_fitness() == 0:
                        search_attribute.append("r")  # change room
                elif three_hours > 0:
                    search_attribute.append("cs")  # change course slot
                elif three_hours < 0:
                    search_attribute.append("c")  # change course
        if parent1.get_prof_handled_course_sched_fitness() == 0 and parent1.get_prof_handled_course_fitness() == 1:
            handled_course = parent1.professor.get_course_code_handle()

            laboratory = False

            course = parent1.get_course().get_code()
            count = handled_course.count(course)

            if course in LAB_COURSE:
                three_hours = 3 * count
                two_hours = 2 * count
                laboratory = True
            else:
                three_hours = 3 * count
                two_hours = 0 * count

            # get all the schedules from monday - saturday
            schedule_course_code = parent1.get_professor().get_schedule_course()

            for sched in schedule_course_code:
                if sched[0] == course:
                    if sched[1] == 3:
                        three_hours -= 3
                    elif sched[1] == 2:
                        two_hours -= 2

                    if three_hours < 0 and two_hours < 0:
                        break

            if laboratory:
                if three_hours == 0 and two_hours > 0:
                    # look for same course with 2 hours time slot
                    search_attribute.append("cs")  # change course slot
                elif three_hours == 0 and two_hours < 0 and parent1.get_course().get_hour() == 2:
                    search_attribute.append("c")  # change course
                elif three_hours > 0 and two_hours == 0:
                    # look for same course with 3 hours time slot
                    search_attribute.append("cs")  # change course slot
                elif three_hours < 0 and two_hours == 0 and parent1.get_course().get_hour() == 3:
                    search_attribute.append("c")  # change course
                elif three_hours < 0 and two_hours < 0:
                    search_attribute.append("c")  # change course
                elif three_hours == 0 and two_hours == 0:
                    if parent1.get_block_enrolled_course_fitness() == 0:
                        search_attribute.append("b")  # change prof
                    if parent1.get_room_suitability_fitness() == 0:
                        search_attribute.append("r")  # change room
            else:
                if three_hours == 0:
                    if parent1.get_block_enrolled_course_fitness() == 0:
                        search_attribute.append("b")  # change prof
                    if parent1.get_room_suitability_fitness() == 0:
                        search_attribute.append("r")  # change room
                elif three_hours > 0:
                    search_attribute.append("cs")  # change course slot
                elif three_hours < 0:
                    search_attribute.append("c")  # change course

        for i in check_attribute:
            if i not in search_attribute:
                break
        else:
            remove_duplicate = list(set(search_attribute))
            search_attribute = remove_duplicate
            random_num = np.random.randint(len(check_attribute))
            search_attribute.remove(check_attribute[random_num])

            if parent1.get_working_hours_fitness() == 0:
                search_attribute.append('d')

        return search_attribute


    def parent_two(self, items, attribute, parent1):

        perfect_gene_found = None
        perfect_gene_score = 0
        prev_time_slot = parent1.get_class_slot_one(2)
        prev_day = parent1.get_class_slot_one(1)

        orig_block_sched = copy.copy(parent1.get_block().get_schedule().get_schedule(prev_day))
        orig_prof_sched = copy.copy(parent1.get_professor().get_schedule().get_schedule(prev_day))
        orig_prof_course_sched = copy.copy(parent1.get_professor().get_schedule_course())
        orig_block_course_sched = copy.copy(parent1.get_block().get_schedule_course())

        for i in range(18):
            random_num = np.random.randint(self.chromosome_size)

            item = items[random_num]

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

                if parent1_copy.get_course().get_code() != item.get_course().get_code():

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
                # return parent2.append([item, 1])
                return item
            elif fitness_score > perfect_gene_score:
                perfect_gene_found = item
                perfect_gene_score = fitness_score

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
                        gene_fitness.block_enrolled_course() +
                        gene_fitness.block_handled_course_schedule() +
                        gene_fitness.prof_handled_course_schedule()
                ) / 14
        )

        return fitness_score

    def uniform_crossover(self, parent1, parent2):
        attribute = GeneticAlgorithm.check_attribute(parent1)

        parent1_copy = copy.copy(parent1)
        parent2_copy = copy.copy(parent2)

        prev_day = parent1_copy.get_class_slot_one(1)

        remove_prof_sched = False
        remove_block_sched = False

        if 'cs' in attribute or 'r' in attribute or 'd' in attribute:

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

            remove_prof_sched = True
            remove_block_sched = True

        if 'p' in attribute:

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
            child = self.mutation(parent1)
            self.add_to_new_chromosome(child)
        elif parent1.get_fitness_score() < parent2.get_fitness_score():
            parent2_copy = copy.copy(parent2)
            self.curr_chromosome.pop_gene(self.curr_chromosome.get_position(parent2))
            self.curr_chromosome.add_gene(parent1)
            child = self.mutation(parent2_copy)
            self.add_to_new_chromosome(child)
        else:
            child = self.mutation(parent1)
            self.add_to_new_chromosome(child)

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

            prev_time_slot_rg = rg_copy.get_class_slot_one(2)
            prev_day_rg = rg_copy.get_class_slot_one(1)

            if 'cs' in attribute or 'r' in attribute or 'd' in attribute:

                if len(self.class_slot_population) > 0:

                    random_class_slot = self.class_slot_population.pop(
                        random.randrange(len(self.class_slot_population)))

                    child_class_slot = child.get_class_slot()
                    self.class_slot_population.append(child_class_slot)

                    child.set_class_slot(random_class_slot)

                    new_time_slot = child.get_class_slot_one(2)
                    new_day = child.get_class_slot_one(1)

                    child_copy.get_class_slot_one(0).get_schedule().remove_sched(prev_day, prev_time_slot)
                    child.get_class_slot_one(0).set_schedule(new_day, new_time_slot)

                else:
                    child_class_slot = child.get_class_slot()
                    child.set_class_slot(random_gene.get_class_slot())
                    random_gene.set_class_slot(child_class_slot)
                    random_gene_picked = True

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)

                if random_gene_picked:
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

                remove_block_sched = True
                remove_prof_sched = True

            if 'p' in attribute:

                child_course = child.get_course().get_code()
                child_course_hour = child.get_course().get_hour()

                rg_course = random_gene.get_course().get_code()
                rg_course_hour = random_gene.get_course().get_hour()

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)
                rg_new_time_slot = random_gene.get_class_slot_one(2)
                rg_new_day = random_gene.get_class_slot_one(1)

                child_prof = child.get_professor()

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

                # magpapalit ang course handle ng prof ng p1 at p2
                child_copy.get_professor().remove_schedule_course(child_course, child_course_hour)
                child.get_professor().set_schedule_course(child_course, child_course_hour)
                rg_copy.get_professor().remove_schedule_course(rg_course, rg_course_hour)
                random_gene.get_professor().set_schedule_course(rg_course, rg_course_hour)

            if 'b' in attribute:

                child_course = child.get_course().get_code()
                child_course_hour = child.get_course().get_hour()

                rg_course = random_gene.get_course().get_code()
                rg_course_hour = random_gene.get_course().get_hour()

                new_time_slot = child.get_class_slot_one(2)
                new_day = child.get_class_slot_one(1)
                rg_new_time_slot = random_gene.get_class_slot_one(2)
                rg_new_day = random_gene.get_class_slot_one(1)

                child_block = child.get_block()
                random_gene_block = random_gene.get_block()

                child.set_block(random_gene_block)
                random_gene.set_block(child_block)

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


                # magpapalit ang course handle ng prof ng p1 at p2
                child_copy.get_block().remove_schedule_course(child_course, child_course_hour)
                child.get_block().set_schedule_course(child_course, child_course_hour)
                rg_copy.get_block().remove_schedule_course(rg_course, rg_course_hour)
                random_gene.get_block().set_schedule_course(rg_course, rg_course_hour)

            if 'c' in attribute:
                child_course = child.get_course().get_code()
                child_course_hour = child.get_course().get_hour()

                rg_course = random_gene.get_course().get_code()
                rg_course_hour = random_gene.get_course().get_hour()

                course = child.get_course()
                child.set_course(random_gene.get_course())
                random_gene.set_course(course)

                # magpapalit ang course handle ng block ng p1 at p2
                child.get_block().remove_schedule_course(child_course, child_course_hour)
                child.get_block().set_schedule_course(rg_course, rg_course_hour)
                random_gene.get_block().remove_schedule_course(rg_course, rg_course_hour)
                random_gene.get_block().set_schedule_course(child_course, child_course_hour)

                child.get_professor().remove_schedule_course(child_course, child_course_hour)
                child.get_professor().set_schedule_course(rg_course, rg_course_hour)
                random_gene.get_professor().remove_schedule_course(rg_course, rg_course_hour)
                random_gene.get_professor().set_schedule_course(child_course, child_course_hour)

        else:
            if 'cs' in attribute or 'r' in attribute or 'd' in attribute:

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
                    self.gene.set_professor_work_load_fitness(sublist[1])
                    return sublist[1]

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

    def prof_handled_course_schedule(self):
        prof_code = self.gene.professor.get_prof_id()
        handled_course = self.gene.professor.get_course_code_handle()

        if MEMOIZATION:
            for sublist in PHCS_CACHE:
                if prof_code == sublist[0]:
                    self.gene.set_prof_handled_course_sched_fitness(sublist[1])
                    return sublist[1]

        counted_elements = Counter(handled_course)

        result = []
        for course, count in counted_elements.items():

            if course in LAB_COURSE:
                three_hours = 3 * count
                two_hours = 2 * count
            else:
                three_hours = 3 * count
                two_hours = 0

            result.append([course, three_hours, two_hours])

        # get all the schedules from monday - saturday
        schedule_course_code = self.gene.professor.get_schedule_course()

        for sched in schedule_course_code:
            for course in result:
                if sched[0] == course[0]:
                    if sched[1] == 3:
                        course[1] -= 3
                    elif sched[1] == 2:
                        course[2] -= 2

        for course in result:
            if course[1] != 0 and course[2] != 0:
                self.gene.set_prof_handled_course_sched_fitness(0)
                self.gene.set_prof_handled_course_sched_fitness(0)
                if MEMOIZATION:
                    PHCS_CACHE.append([prof_code, 0])
                return 0

        if MEMOIZATION:
            PHCS_CACHE.append([prof_code, 1])

        self.gene.set_prof_handled_course_sched_fitness(1)

        return 1

    def block_handled_course_schedule(self):
        block_code = self.gene.block.get_block_code()
        handled_course = self.gene.block.get_courses()

        if MEMOIZATION:
            for sublist in BHCS_CACHE:
                if block_code == sublist[0]:
                    self.gene.set_block_handled_course_sched_fitness(sublist[1])
                    return sublist[1]

        counted_elements = Counter(handled_course)

        result = []
        for course, count in counted_elements.items():

            if course in LAB_COURSE:
                three_hours = 3 * count
                two_hours = 2 * count
            else:
                three_hours = 3 * count
                two_hours = 0

            result.append([course, three_hours, two_hours])

        # get all the schedules from monday - saturday
        schedule_course_code = self.gene.professor.get_schedule_course()

        for sched in schedule_course_code:
            for course in result:
                if sched[0] == course[0]:
                    if sched[1] == 3:
                        course[1] -= 3
                    elif sched[1] == 2:
                        course[2] -= 2

        for course in result:
            if course[1] != 0 and course[2] != 0:
                self.gene.set_prof_handled_course_sched_fitness(0)
                if MEMOIZATION:
                    BHCS_CACHE.append([block_code, 0])
                return 0

        self.gene.set_block_handled_course_sched_fitness(1)
        if MEMOIZATION:
            BHCS_CACHE.append([block_code, 1])

        return 1


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
        self.schedule_course = []
        self.schedule = Schedule()

    def get_prof_id(self):
        return self.prof_id

    def get_course_code_handle(self):
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

    sl1 = ClassRoom("SL1", 50, "Laboratory")
    sl2 = ClassRoom("SL2", 50, "Laboratory")
    sl3 = ClassRoom("SL3", 50, "Laboratory")
    cisco = ClassRoom("CISCO", 50, "Laboratory")
    r1 = ClassRoom("101", 50, "Regular")
    r2 = ClassRoom("102", 50, "Regular")
    r3 = ClassRoom("103", 50, "Regular")
    r4 = ClassRoom("104", 50, "Regular")

    Data.add_room(sl1)
    Data.add_room(sl2)
    Data.add_room(sl3)
    Data.add_room(cisco)
    Data.add_room(r1)
    Data.add_room(r2)
    Data.add_room(r3)
    Data.add_room(r4)

    # year = 1 if first year else 0
    Data.add_lab_course(LabCourse("CS 121", 3, ["Sl1", "SL2", "SL3"]))
    Data.add_lab_course(LabCourse("CS 222", 3, ["Sl1", "SL2", "SL3", "CISCO"]))
    Data.add_lab_course(LabCourse("CS 221", 3, ["Sl1", "SL2"]))
    Data.add_lab_course(LabCourse("CS 322", 3, ["Sl1", "SL2", "SL3"]))

    Data.add_lab_course(LabCourse("CS 121", 2, ["101", "102", "103", "104"]))
    Data.add_lab_course(LabCourse("CS 222", 2, ["101", "102", "103", "104"]))
    Data.add_lab_course(LabCourse("CS 221", 2, ["101", "102", "103", "104"]))
    Data.add_lab_course(LabCourse("CS 322", 2, ["101", "102", "103", "104"]))

    Data.add_lec_course(LectureCourse("MATH 111", 3, ["101", "102", "103", "104"]))
    Data.add_lec_course(LectureCourse("GEd 105", 3, ["101", "102", "103", "104"]))
    Data.add_lec_course(LectureCourse("GEd 107", 3, ["101", "102", "103", "104"]))
    Data.add_lec_course(LectureCourse("CS 323", 3, ["101", "102", "103", "104"]))

    block101 = Block("101", 45, ["CS 121", "CS 222", "MATH 111"], 1)
    block102 = Block("102", 45, ["CS 121", "CS 221", "MATH 111"], 1)
    block201 = Block("201", 45, ["CS 322", "CS 323", "GEd 107"], 0)
    block202 = Block("202", 45, ["GEd 105", "CS 222", "CS 121"], 0)
    block203 = Block("203", 45, ["GEd 105", "CS 222", "CS 121"], 0)
    block301 = Block("301", 45, ["CS 322", "CS 323", "GEd 107"], 0)

    Data.add_block(block101)
    Data.add_block(block102)
    Data.add_block(block201)
    Data.add_block(block202)
    Data.add_block(block203)
    Data.add_block(block301)

    p1 = Professor("P1", ["CS 121", "CS 222", "MATH 111"])
    p2 = Professor("P2", ["CS 121", "CS 221", "CS 323"])
    p3 = Professor("P3", ["MATH 111", "CS 322", "GEd 107"])
    p4 = Professor("P4", ["GEd 105", "CS 322", "CS 222", "CS 323"])
    p5 = Professor("P5", ["CS 121", "GEd 107"])
    p6 = Professor("P6", ["CS 121", "CS 222", "GEd 105"])

    Data.add_prof(p1)
    Data.add_prof(p2)
    Data.add_prof(p3)
    Data.add_prof(p4)
    Data.add_prof(p5)
    Data.add_prof(p6)

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
            prof = [p1, p2, p3, p4, p5, p6]
            room = [sl1, sl2, sl3, cisco, r1, r2, r3, r4]
            block = [block101, block102, block201, block202, block203, block301]

            for i in range(len(prof)):
                print("===========================================")
                print("Professor ", i+1)
                print("Monday : ", prof[i].get_schedule().get_m())
                print("Tuesday : ", prof[i].get_schedule().get_t())
                print("Wednesday : ", prof[i].get_schedule().get_w())
                print("Thursday : ", prof[i].get_schedule().get_th())
                print("Friday : ", prof[i].get_schedule().get_f())
                print("Saturday : ", prof[i].get_schedule().get_s())
                print("===========================================")

            for i in range(len(block)):
                print("===========================================")
                print("Block ", i + 1)
                print("Monday : ", block[i].get_schedule().get_m())
                print("Tuesday : ", block[i].get_schedule().get_t())
                print("Wednesday : ", block[i].get_schedule().get_w())
                print("Thursday : ", block[i].get_schedule().get_th())
                print("Friday : ", block[i].get_schedule().get_f())
                print("Saturday : ", block[i].get_schedule().get_s())
                print("===========================================")

            for i in range(len(room)):
                print("===========================================")
                print("Room ", i + 1)
                print("Monday : ", room[i].get_schedule().get_m())
                print("Tuesday : ", room[i].get_schedule().get_t())
                print("Wednesday : ", room[i].get_schedule().get_w())
                print("Thursday : ", room[i].get_schedule().get_th())
                print("Friday : ", room[i].get_schedule().get_f())
                print("Saturday : ", room[i].get_schedule().get_s())
                print("===========================================")
            not_perfect_schedule = False
        else:
            while new_chromosome:
                parent1 = genetic_algo.roulette_wheel()

                if parent1 is not None:

                    parent2 = genetic_algo.hill_climbing(parent1)

                    if parent2 is None:
                        genetic_algo.add_to_new_chromosome(parent1)
                    elif parent2 == 0:
                        child = genetic_algo.mutation(parent1)
                        genetic_algo.add_to_new_chromosome(child)
                    else:
                        genetic_algo.uniform_crossover(parent1, parent2)
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