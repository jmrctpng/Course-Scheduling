from typing import List, Tuple, Any
import random


class Data:
    professor = []
    info_technology = []
    comp_science = []
    lab_room = []
    reg_room = []
    lab_course = []
    reg_course = []

    @staticmethod
    def add_lab_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Data.lab_course.append(course)

    @staticmethod
    def add_room(room, type_of_room):
        if type_of_room == "Regular":
            Data.reg_room.append(room)
        else:
            Data.lab_room.append(room)

    @staticmethod
    def find_prof(prof_id):
        for prof in Data.professor:
            if prof.get_prof_id() == prof_id:
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
    def find_block(block_code):
        blocks = Data.info_technology + Data.comp_science

        for block in blocks:
            if block.get_block_code() == block_code:
                return block

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


class GeneticAlgorithm:

    def __init__(self):
        self.chromosome = Chromosome()

    def encode_chromosome(self):
        """
        Randomly select one item from each population list and create a new chromosome
        list containing those four items. The for loop iterates over the length of
        self.course_population and creates one chromosome for each course.
        """

        population_class = Population()
        data = Data()

        for i in range(len(population_class.courses)):
            class_slot = population_class.class_slots.pop(random.randrange(len(population_class.class_slots)))
            course = population_class.courses.pop(random.randrange(len(population_class.courses)))
            block = population_class.blocks.pop(random.randrange(len(population_class.blocks)))
            professor = population_class.professors.pop(random.randrange(len(population_class.professors)))

            gene = Gene(class_slot, course, block, professor)

            # set schedule for professor
            found_prof = data.find_prof(professor)
            if found_prof is not None:
                found_prof.schedule.set_schedule(class_slot[1], class_slot[2])

            # set schedule for classroom
            found_room = data.find_room(class_slot[0].get_code(), class_slot[0].type_of_room())
            if found_room is not None:
                found_room.schedule.set_schedule(class_slot[1], class_slot[2])

            # set schedule for block
            found_block = data.find_block(block.get_block_code())
            if found_block is not None:
                found_block.set_schedule(class_slot[1], class_slot[2])

            self.chromosome.add_gene(gene)

    def calculate_fit(self):
        pass

    def roulette_wheel(self):
        pass


class Chromosome:
    def __init__(self):
        self.genes = []

    def add_gene(self, gene):
        self.genes.append(gene)

    def get_genes(self):
        return self.genes


class Gene:
    def __init__(self, class_slot, course, block, professor):
        self.class_sot = class_slot
        self.course = course
        self.block = block
        self.professor = professor

    def get_class_slot(self):
        return self.class_sot

    def get_class_slot_one(self, index):
        return self.class_sot[index]

    def get_class_slot_two(self, index1, index2):
        return self.class_sot[index1][index2]

    def get_course(self):
        return self.course

    def get_block(self):
        return self.block

    def get_professor(self):
        return self.professor

    def get_fitness_score(self):
        pass


class GeneFitness:

    def __init__(self, gene):
        self.gene = gene

    def classroom_capacity(self):
        if self.gene.block.get_enrolled_students() <= self.gene.get_class_slot_one(0).get_capacity():
            return 1
        else:
            return 0

    def professor_work_load(self):
        day = self.gene.get_class_slot_one(1)
        schedule = self.gene.professor.get_schedule()

        day_sched = getattr(schedule, f'get_{day.lower()}')()

        result = has_overlap(self.gene.get_class_slot_two(2, 0), day_sched)

        return -1 if result else 0

    def schedule_availability(self):
        pass

    def class_room_availability(self):
        room = self.gene.get_class_slot_one(0,0)
        schedule = room.get_schedule()

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



    def room_suitability(self):
        room = self.gene.get_class_slot_one(0, 0)
        available_rooms = self.gene.course.get_available_rooms()

        if room in available_rooms:
            return 1
        else:
            return 0

    def first_year(self):
        if self.gene.course.get_year() == 1:
            if self.gene.get_class_slot_one(1) != 'S':
                return 1
            else:
                return 0

    def has_lunch_break(self):

        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        # Get all the sublists that contain a number between 10 and 14
        schedules = [schedule for schedule in professor_schedule if any(10 <= int(time) <= 14 for time in schedule)]

        # Create a new list with a sublist that contains only two elements
        relevant_sublists = []
        for schedule in schedules:
            for i in range(len(schedule) - 1):
                sublist = [int(schedule[i]), int(schedule[i + 1])]
                if 10 <= sublist[0] <= 14 and 10 <= sublist[1] <= 14:
                    relevant_sublists.append(sublist)

        # Sort the array
        sorted_sublists = sorted(relevant_sublists)

        # Check if the list of relevant sublists is a subset of the sorted list
        return set([tuple(sublist) for sublist in [[10, 11], [11, 12], [12, 13], [13, 14]]]).issubset(
            [tuple(sublist) for sublist in sorted_sublists])

    def maximum_working_hours(self):
        prof_sched_day = self.gene.get_class_slot_one(1)
        professor_schedule = self.gene.professor.schedule.get_schedule(prof_sched_day)

        total = 0
        for sublist in professor_schedule:
            length = len(sublist) - 1
            total += length

        if total > 8:
            total = total - 8

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
        return self.capacity()

    def get_schedule(self):
        return self.schedule

    def get_type_of_room(self):
        return self.type_of_room


class Course:

    def __init__(self, code, lab_hours, lecture_hours, available_rooms, year):
        self.code = code
        self.lab_hours = lab_hours
        self.lecture_hours = lecture_hours
        self.available_rooms = available_rooms
        self.year = year

    def get_code(self):
        return self.code()

    def get_lab_hours(self):
        return self.lab_hours()

    def get_lecture_hours(self):
        return self.lecture_hours()

    def get_available_rooms(self):
        return self.available_rooms

    def get_year(self):
        return self.year


class Professor:
    schedule = Schedule()

    def __init__(self, prof_id, course_code_handle):
        self.prof_id = prof_id
        self.course_code_handle = course_code_handle

    def get_prof_id(self):
        return self.prof_id()

    def get_course_code_handle(self):
        return self.course_code_handle()

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

    def initialize_population(self, classSlot, professor, course, block):
        pass

    def generate_course_population(self):
        block = Data()
        available_course = block.comp_science + block.info_technology

        for b in available_course:
            for course in b.get_courses():
                self.courses.append(course)

    def generate_block_population(self):
        block = Data()
        available_blocks = block.comp_science + block.info_technology

        for b in available_blocks:
            for i in range(len(b.get_courses())):
                self.blocks.append(b.get_block_code())

    def generate_prof_population(self):
        professor = Data()

        for prof in professor.professor:
            for i in range(len(prof.get_course_code_handle())):
                self.professors.append(prof.get_prof_id())

    def generate_class_slot_population(self, meetingTime, day):
        class_slot = labRoomMapping(meetingTime, day) + regRoomMapping(day)
        self.class_slots.append(class_slot)

    def get_population(self):
        return self.population


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


if __name__ == '__main__':
    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F", "S"]

    Data.add_lab_room("SL1", 50)
    Data.add_lab_room("SL2", 50)
    Data.add_lab_room("SL3", 50)
    Data.add_lab_room("CISCO", 50)
    Data.add_lab_room("ML", 50)
    Data.add_lab_room("ITL", 50)
    Data.add_reg_room("101", 50)
    Data.add_reg_room("102", 50)
    Data.add_reg_room("103", 50)
    Data.add_reg_room("104", 50)
    Data.add_reg_room("201", 50)
    Data.add_reg_room("202", 50)
    Data.add_reg_room("203", 50)
    Data.add_reg_room("EDL", 50)
    Data.add_reg_room("401", 50)
    Data.add_reg_room("402", 50)
    Data.add_reg_room("403", 50)

    Data.add_lab_course("CS 121", 2, 3, ["Sl1", "SL2", "SL3"])
    Data.add_lab_course("IT 221", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("CS 222", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    Data.add_lab_course("CS 221", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    Data.add_lab_course("CS 322", 2, 3, ["SL1", "SL2", "SL3"])
    Data.add_lab_course("CS 324", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("CSE 401", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("CS 422", 2, 3, ["ML"])
    Data.add_lab_course("IT 323", 2, 3, ["ML"])
    Data.add_lab_course("CSE 403", 2, 3, ["ML"])
    Data.add_lab_course("CS 111", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    Data.add_lab_course("IT 221", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("IT 223", 2, 3, ["CISCO"])
    Data.add_lab_course("IT 222", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("NTT 403", 2, 3, ["CISCO"])
    Data.add_lab_course("NTT 404", 2, 3, ["CISCO"])
    Data.add_lab_course("IT 322", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
    Data.add_lab_course("IT 323", 2, 3,
                        ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                         "ML", "ITL", "EDL", "401", "402", "403"])
