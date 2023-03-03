from typing import List, Tuple, Any
import random


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
        blocks = self.schedule[day]
        blocks.append(time)
        self.schedule[day] = blocks

    def get_schedule(self, day):
        return self.schedule[day]

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

    def encode_chromosome(self, population_class):
        """
        Randomly select one item from each population list and create a new chromosome
        list containing those four items. The for loop iterates over the length of
        self.course_population and creates one chromosome for each course.
        """

        data = Data()

        for i in range(len(population_class.get_courses())):
            class_slot = population_class.get_class_slots().pop(random.randrange(len(population_class.class_slots)))
            course = population_class.get_courses().pop(random.randrange(len(population_class.courses)))
            block = population_class.get_blocks().pop(random.randrange(len(population_class.blocks)))
            professor = population_class.get_professors().pop(random.randrange(len(population_class.professors)))

            gene = Gene(class_slot, course, block, professor)

            # set schedule for professor
            found_prof = data.find_prof(professor)
            if found_prof is not None:
                found_prof.set_schedule(class_slot[1], class_slot[2])

            # set schedule for classroom
            found_room = data.find_room(class_slot[0].get_code(), class_slot[0].get_type_of_room())
            if found_room is not None:
                found_room.set_schedule(class_slot[1], class_slot[2])

            # set schedule for block
            found_block = data.find_block(block)
            if found_block is not None:
                found_block.set_schedule(class_slot[1], class_slot[2])

            self.chromosome.add_gene(gene)

        GeneticAlgoPopulation.chromosome.append(self.chromosome)

    @staticmethod
    def calculate_fit():
        # get the current chromosome
        if len(GeneticAlgoPopulation.chromosome) == 1:
            chromosome = GeneticAlgoPopulation.chromosome[0]
        else:
            chromosome = GeneticAlgoPopulation.chromosome[-1]

        # calculate gene fitness score
        for gene in chromosome.get_genes():
            gene_fitness = GeneFitness(gene)
            n = 10
            first_year = 0

            if gene_fitness.first_year() is None:
                n = 9
            else:
                first_year = 1


            fitness_score = (
                (
                    gene_fitness.classroom_capacity() +
                    gene_fitness.professor_work_load() +
                    gene_fitness.schedule_availability() +
                    gene_fitness.room_suitability() +
                    first_year +
                    gene_fitness.block_lunch() +
                    gene_fitness.professor_lunch() +
                    gene_fitness.maximum_working_hours()
                 ) / n
            )

            gene.add_fitness_score(fitness_score)
            print("Class slot: ", gene.class_slot[0].get_code(), " ", gene.class_slot[0].get_capacity(), " ", gene.class_slot[0].get_type_of_room(), " ", gene.class_slot[1], " ", gene.class_slot[2])
            print("Course: ", gene.course)
            print("Professor: ", gene.professor.get_prof_id())
            print("Block: ", gene.block.get_block_code())
            print("Score :", gene.fitness_score)
            print("CR :", gene.get_classroom_capacity_fitness())
            print("PW :", gene.get_professor_work_load_fitness())
            print("SA :", gene.get_room_availability_fitness())
            print("PW :", gene.get_professor_availability_fitness())
            print("BA :", gene.get_block_availability_fitness())
            print("RS :", gene.get_room_suitability_fitness())
            print("FF :", first_year)
            print("PL :", gene.get_prof_lunch_break_fitness())
            print("BL :", gene.get_block_lunch_break_fitness())
            print("WH :", gene.get_working_hours_fitness())
            print("-----------------------------")

    def roulette_wheel(self):
        pass


class Chromosome:
    def __init__(self):
        self.genes = []

    def add_gene(self, gene):
        self.genes.append(gene)

    def get_genes(self):
        return self.genes

    def get_specific_gene(self, index):
        return self.genes[index]


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
        self.professor_availability_fitness = 0
        self.block_availability_fitness = 0
        self.room_suitability_fitness = 0
        self.first_year_fitness = 0
        self.prof_lunch_break_fitness = 0
        self.block_lunch_break_fitness = 0
        self.working_hours_fitness = 0

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

    def get_professor_availability_fitness(self):
        return self.professor_availability_fitness

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

    def set_classroom_capacity_fitness(self, value):
        self.classroom_capacity_fitness = value

    def set_professor_work_load_fitness(self, value):
        self.professor_work_load_fitness = value

    def set_room_availability_fitness(self, value):
        self.room_availability_fitness = value

    def set_professor_availability_fitness(self, value):
        self.professor_availability_fitness = value

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

        #day_sched = getattr(schedule, f'get_{day.lower()}')()

        result = has_overlap(self.gene.get_class_slot_two(2, 0), schedule)

        if result:
            self.gene.set_professor_availability_fitness(0)
            return -1
        else:
            self.gene.set_professor_availability_fitness(1)
            return 1

    def schedule_availability(self):

        day = self.gene.get_class_slot_one(1)
        room_schedule = self.gene.class_slot[0].get_schedule().get_schedule(day)
        block_schedule = self.gene.block.get_schedule().get_schedule(day)
        prof_schedule = self.gene.professor.get_schedule().get_schedule(day)

        room_availability = True if available_hours_in_schedule(room_schedule) > 2 else False
        block_availability = True if available_hours_in_schedule(block_schedule) > 2 else False
        prof_availability = True if available_hours_in_schedule(prof_schedule) > 2 else False

        self.gene.set_room_availability_fitness(1) if room_availability else 0
        self.gene.set_professor_availability_fitness(1) if prof_availability else 0
        self.gene.set_block_availability_fitness(1) if block_availability else 0

        if room_availability and block_availability and prof_availability:
            return 3
        if room_availability + block_availability + prof_availability == 2:
            return 2
        if any([room_availability, block_availability, prof_availability]):
            return 1
        return 0

    def room_suitability(self):
        room = self.gene.get_class_slot_one(0)
        course = Data.find_course(self.gene.course)
        available_rooms = course.get_available_rooms()

        if room.get_code() in available_rooms:
            self.gene.set_room_suitability_fitness(1)
            return 1
        else:
            self.gene.set_room_suitability_fitness(0)
            return 0

    def first_year(self):

        course = Data.find_course(self.gene.course)

        if course.get_year() == 1:
            if self.gene.get_class_slot_one(1) != 'S':
                self.gene.set_first_year_fitness(1)
                return 1
            else:
                self.gene.set_first_year_fitness(0)
                return 0
        else:
            self.gene.set_first_year_fitness(None)
            return None

    def professor_lunch(self):
        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        # if true then no lunch
        if has_no_lunch_break(professor_schedule):
            self.gene.set_prof_lunch_break_fitness(0)
            return 0
        else:
            self.gene.set_prof_lunch_break_fitness(1)
            return 1

    def block_lunch(self):
        block_sched_day = self.gene.get_class_slot_one(1)
        block_schedule = self.gene.block.schedule.get_schedule(block_sched_day)

        # if true then no lunch
        if has_no_lunch_break(block_schedule):
            self.gene.set_block_lunch_break_fitness(0)
            return 0
        else:
            self.gene.set_block_lunch_break_fitness(1)
            return 1

    def maximum_working_hours(self):
        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        total = 0
        for sublist in professor_schedule:
            length = len(sublist) - 1
            total += length

        if total > 8:
            total = 8 - total
            self.gene.set_working_hours_fitness(0)
        else:
            self.gene.set_working_hours_fitness(1)
            total = 0

        return total


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

    def __init__(self, code, lab_hours, available_rooms, year):
        self.code = code
        self.lab_hours = lab_hours
        self.available_rooms = available_rooms
        self.year = year

    def get_code(self):
        return self.code

    def get_lab_hours(self):
        return self.lab_hours()

    def get_available_rooms(self):
        return self.available_rooms

    def get_year(self):
        return self.year


class LectureCourse:

    def __init__(self, code, lecture_hours, available_rooms, year):
        self.code = code
        self.lecture_hours = lecture_hours
        self.available_rooms = available_rooms
        self.year = year

    def get_code(self):
        return self.code()

    def get_lecture_hours(self):
        return self.lecture_hours()

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
    schedule = Schedule()

    def __init__(self, block_code, enrolled_students, courses):
        self.block_code = block_code
        self.enrolled_students = enrolled_students
        self.courses = courses

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
        # available_course = data.lab_course + data.lec_course
        blocks = data.block

        for block in blocks:
            for course in block.get_courses():
                self.courses.append(course)
                self.blocks.append(block)

    def generate_prof_population(self):
        data = Data()

        for prof in data.professor:
            for course in prof.get_course_code_handle():
                for data_lec_course in Data.lec_course:
                    if course == data_lec_course.get_code():
                        self.professors.append(prof)
                        break
                for data_lab_course in Data.lab_course:
                    if course == data_lab_course.get_code():
                        self.professors.append(prof)
                        break


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


def has_overlap(time_slot, schedule):
    for slot in time_slot:
        for i, time in enumerate(schedule):
            if slot[0] in time and i < len(schedule) - len(slot) + 1:
                if all(slot[j] == schedule[i + j] for j in range(1, len(slot))):
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
    Data.add_lab_course(LabCourse("CS 121", 3, ["Sl1", "SL2", "SL3"], 1))
    Data.add_lab_course(LabCourse("CS 222", 3, ["Sl1", "SL2", "SL3", "CISCO"], 0))
    Data.add_lab_course(LabCourse("CS 221", 3, ["Sl1", "SL2"], 0))
    Data.add_lab_course(LabCourse("CS 322", 3, ["Sl1", "SL2", "SL3"], 0))

    Data.add_lec_course(LabCourse("CS 121", 2, ["Sl1", "SL2", "SL3"], 1))
    Data.add_lec_course(LabCourse("CS 222", 2, ["Sl1", "SL2", "SL3", "CISCO"], 0))
    Data.add_lec_course(LabCourse("CS 221", 2, ["Sl1", "SL2"], 0))
    Data.add_lec_course(LabCourse("CS 322", 2, ["Sl1", "SL2", "SL3"], 0))

    Data.add_lec_course(LabCourse("MATH 111", 3, ["101", "102", "103", "104"], 0))
    Data.add_lec_course(LabCourse("GEd 105", 3, ["101", "102", "103", "104"], 0))
    Data.add_lec_course(LabCourse("GEd 107", 3, ["101", "102", "103", "104"], 0))
    Data.add_lec_course(LabCourse("CS 323", 3, ["101", "102", "103", "104"], 0))

    Data.add_block(Block("101", 45, ["CS 121", "CS 222", "MATH 111"]))
    Data.add_block(Block("102", 45, ["CS 121", "CS 221", "MATH 111"]))
    Data.add_block(Block("201", 45, ["CS 322", "CS 323", "GEd 107"]))
    Data.add_block(Block("202", 45, ["GEd 105", "CS 222", "CS 121"]))
    Data.add_block(Block("203", 45, ["GEd 105", "CS 222", "CS 121"]))
    Data.add_block(Block("301", 45, ["CS 322", "CS 323", "GEd 107"]))

    Data.add_prof(Professor("P1", ["CS 121", "CS 222", "MATH 111"]))
    Data.add_prof(Professor("P2", ["CS 121", "CS 221", "CS 323"]))
    Data.add_prof(Professor("P3", ["MATH 111", "CS 322", "GEf 107"]))
    Data.add_prof(Professor("P4", ["GEd 105", "CS 322", "CS 222", "CS 323"]))
    Data.add_prof(Professor("P5", ["CS 121", "GEd 107"]))

    population = Population()
    population.initialize_population(meetingTime, day)

    genetic_algo = GeneticAlgorithm()
    genetic_algo.encode_chromosome(population)
    genetic_algo.calculate_fit()

