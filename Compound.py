from dataCollection import periodicTable

class Compound:
    def __init__(self, equation):
        self.charge = 0  # Default charge is 0 for neutral compounds
        self.equation = equation
        self.elements = self.parseElements(equation)

        self.elements.setdefault("e-", 0)
        self.charge = - self.elements["e-"]

        del self.elements["e-"]

    def __str__(self):
        return f"{self.equation} - charge {self.charge}"

    def __repr__(self):
        return f"Compound(equation={self.equation}, charge={self.charge})"
    
    def parseElements(self, equation):
        if len(equation) == 0: return {}
        if equation[0] == "_": return self.parseElements(equation[1:])
        if equation[0] == "+": 
            if len(equation) == 1: return {"e-" : 1}
            return {"e-" : -1 * int(equation[1:])}
        if equation[0] == "-": 
            if len(equation) == 1: return {"e-" : -1}
            return {"e-" : int(equation[1:])}
        if checkLast(equation): return {equation : 1}
        if equation[-1].isnumeric() and checkLast(equation[:-1]): return {equation[:-1] : int(equation[-1])}
        if equation[-2:].isnumeric() and checkLast(equation[:-2]): return {equation[:-2] : int(equation[-2:])}
        if equation[:5] == "(NH4)":
            factor = int(equation[5])
            curr_batch = {"N" : factor, "H" : 4 * factor}
            curr_batch.update(self.parseElements(equation[6:]))
        if equation[0] == "(":
            factor = int(equation[-1])
            inside = self.parseElements(equation[1:-2])
            for key in inside: inside[key] *= factor
            return inside

        index = int(equation[1].islower())
        num = equation[index+1]
        num = 1 if not num.isnumeric() else int(num)

        ret = {equation[0:index + 1] : num}
        next_batch = self.parseElements(equation[index + 2 - int(num == 1):])
        to_delete = []
        for key in next_batch:
            if key in ret: 
                to_delete.append(key)
                ret[key] += next_batch[key]
        for i in to_delete: del next_batch[i]
        ret.update(next_batch)
        return ret
    
    def totalElements(self):
        total = 0
        for value in self.elements.values():
            total += value
        return total - self.elements.get("e-", 0)
    
    def numUniqueElements(self):
        return len(self.elements) - (1 if "e-" in self.elements else 0)
    
    def uniqueElements(self):
        return [key for key in self.elements if key != "e-"]
    
    def getCentralAtom(self):
        options = self.uniqueElements()

        if len(options) == 1: return options[0]

        if "H" in options:
            options.remove("H")

        if len(options) == 0:
            return False

        # Count occurrences of each non-H element
        counts = {el: self.elements[el] for el in options}

        # If more than one non-H element appears multiple times, reject
        multi_center_candidates = [el for el, count in counts.items() if count > 1]
        if len(multi_center_candidates) > 1:
            return False

        sorted_by_en = sorted(options, key=lambda x: periodicTable[x][2])
        lowest_en = periodicTable[sorted_by_en[0]][2]

        # Check how many elements share this electronegativity
        same_en = [el for el in sorted_by_en if periodicTable[el][2] == lowest_en]

        # If multiple elements share lowest electronegativity, it is ambiguous
        if len(same_en) > 1:
            return False

        return sorted_by_en[0]
    
    def isCovalent(self):
        for element in self.uniqueElements():
            if element not in periodicTable:
                return False
            if periodicTable[element][2] < 1.6:  # EN < 1.6 means it is likely metal
                return False

        return True


    

def checkLast(equation : str): return len(equation) == 1 or (equation[-1].islower() and len(equation) == 2)