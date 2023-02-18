from typing import List, Tuple


class ClassRoom:
    lab_rooms = []
    reg_rooms = []

    @staticmethod
    def add_lab_room(code: str, capacity: int):
        ClassRoom.lab_rooms.append((code, capacity))

    @staticmethod
    def add_reg_room(code: str, capacity: int):
        ClassRoom.reg_rooms.append((code, capacity))

    def show_all_lab_rooms(self) -> List[Tuple[str, int]]:
        return self.lab_rooms

    def show_all_reg_rooms(self) -> List[Tuple[str, int]]:
        return self.reg_rooms


class Course:
    lab_courses = []
    reg_courses = []

    @staticmethod
    def add_lab_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Course.lab_courses.append(course)

    @staticmethod
    def add_reg_course(course_code: str, lab_hours: int, lecture_hours: int, available_rooms: List[str]):
        course = [course_code, [lab_hours, lecture_hours], available_rooms]
        Course.reg_courses.append(course)


class Professor:
    professors = []

    @staticmethod
    def add_prof(prof_id: str, course_code_handle: str):
        element = [prof_id, course_code_handle]
        Professor.professors.append(element)


class Program:
    info_technology_list = []
    comp_sci_list = []

    @staticmethod
    def add_info_tech_block(block_code: str, total_students: int, courses: List[str]):
        """
        Adds a block of information technology to the program.
        """
        Program.info_technology_list.append((block_code, total_students, courses))

    @staticmethod
    def add_comp_sci_block(block_code: str, total_students: int, courses: List[str]):
        """
        Adds a block of computer science to the program.
        """
        Program.comp_sci_list.append((block_code, total_students, courses))


class Population:
    population = []
    course_population: List[str] = []
    block_population = []
    professor_population = []
    class_slot_population = []

    @staticmethod
    def createPopulation(classSlot, professor, course, block):
        pass

    @staticmethod
    def encode_chromosome(self):
        pass

    def generate_course_population(self):
        block = Program()
        available_rooms = block.comp_sci_list + block.info_technology_list

        for room in available_rooms:
            for course in room[3]:
                self.course_population.append(course)

    def generate_block_population(self):
        block = Program()
        available_rooms = block.comp_sci_list + block.info_technology_list

        for room in available_rooms:
            for i in range(len(room[3])):
                self.block_population.append(room[i][0])

    def generate_prof_population(self):
        professor = Professor()

        for prof in professor.professors:
            for i in range(len(prof[1])):
                self.professor_population.append(prof[0])

    def generate_class_slot_population(self, meetingTime, day):
        class_slot = labRoomMapping(meetingTime, day) + regRoomMapping(day)
        self.class_slot_population.append(class_slot)

    def getPopulation(self):
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

    labRoom = ClassRoom()
    lab_map = [[lab, day, meetingTime[i:i + 4]]
               for lab in labRoom.lab_rooms
               for day in meetingDays
               for i in range(0, len(meetingTime), 3)
               if len(meetingTime[i:i + 4]) == 4]

    return lab_map


def regRoomMapping(meetingDays: List[str]) -> List[List[str]]:
    reg_room = ClassRoom()
    pattern = [["07", "08", "09"], ["07", "08", "09", "10"], ["09", "10", "11"], ["09", "10", "11", "12"],
               ["10", "11", "12"], ["10", "11", "12", "13"],
               ["11", "12", "13"], ["11", "12", "13", "14"], ["12", "13", "14"], ["12", "13", "14", "15"],
               ["13", "14", "15"], ["13", "14", "15", "16"],
               ["14", "15", "16"], ["14", "15", "16", "17"], ["15", "16", "17"], ["16", "17", "18", "19"],
               ["17", "18", "19"]]

    reg_room_map = [[room, day, meetingTime]
                    for room in reg_room.reg_rooms
                    for day in meetingDays
                    for meetingTime in pattern
                    ]

    return reg_room_map


if __name__ == '__main__':
    meetingTime = ["07", "08", "09", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19"]
    day = ["M", "T", "W", "Th", "F", "S"]

    available_room = ClassRoom()
    available_room.add_lab_room("SL1", 50)
    available_room.add_lab_room("SL2", 50)
    available_room.add_lab_room("SL3", 50)
    available_room.add_lab_room("CISCO", 50)
    available_room.add_lab_room("ML", 50)
    available_room.add_lab_room("ITL", 50)
    available_room.add_reg_room("101", 50)
    available_room.add_reg_room("102", 50)
    available_room.add_reg_room("103", 50)
    available_room.add_reg_room("104", 50)
    available_room.add_reg_room("201", 50)
    available_room.add_reg_room("202", 50)
    available_room.add_reg_room("203", 50)
    available_room.add_reg_room("EDL", 50)
    available_room.add_reg_room("401", 50)
    available_room.add_reg_room("402", 50)
    available_room.add_reg_room("403", 50)

    available_course = Course()
    available_course.add_lab_course("CS 121", 2, 3, ["Sl1", "SL2", "SL3"])
    available_course.add_lab_course("IT 221", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("CS 222", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    available_course.add_lab_course("CS 221", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    available_course.add_lab_course("CS 322", 2, 3, ["SL1", "SL2", "SL3"])
    available_course.add_lab_course("CS 324", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("CSE 401", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("CS 422", 2, 3, ["ML"])
    available_course.add_lab_course("IT 323", 2, 3, ["ML"])
    available_course.add_lab_course("CSE 403", 2, 3, ["ML"])
    available_course.add_lab_course("CS 111", 2, 3, ["SL1", "SL2", "SL3", "ML"])
    available_course.add_lab_course("IT 221", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("IT 223", 2, 3, ["CISCO"])
    available_course.add_lab_course("IT 222", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("NTT 403", 2, 3, ["CISCO"])
    available_course.add_lab_course("NTT 404", 2, 3, ["CISCO"])
    available_course.add_lab_course("IT 322", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
    available_course.add_lab_course("IT 323", 2, 3,
                                    ["101", "102", "103", "104", "201", "202", "203", "SL1", "SL2", "SL3", "CISCO",
                                     "ML", "ITL", "EDL", "401", "402", "403"])
