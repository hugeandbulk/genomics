
def DFS_dist_from_node(query_node, parents):
    """Return dictionary containing distances of parent GO nodes from the query"""
    result = {}
    stack = []
    stack.append( (query_node, 0) )
    while len(stack) > 0:
        print("stack=", stack)
        node, dist = stack.pop()
        result[node] = dist
        if node in parents:
            for parent in parents[node]:
                # Get the first member of each tuple, see 
                # http://stackoverflow.com/questions/12142133/how-to-get-first-element-in-a-list-of-tuples
                stack_members = [x[0] for x in stack]
                if parent not in stack_members:
                    stack.append( (parent, dist+1) )
    return result
    
def BFS_dist_from_node(query_node, parents):
    """Return dictionary containing minimum distances of parent GO nodes from the query"""
    result = {}
    queue = []
    queue.append( (query_node, 0) )
    while queue:
        print("queue=", queue)
        node, dist = queue.pop(0)
        result[node] = dist
        if node in parents: # If the node *has* parents
            for parent in parents[node]:
                # Get the first member of each tuple, see 
                # http://stackoverflow.com/questions/12142133/how-to-get-first-element-in-a-list-of-tuples
                queue_members = [x[0] for x in queue]
                if parent not in result and parent not in queue_members: # Don't visit a second time
                    queue.append( (parent, dist+1) )
    return result
    
if __name__ == "__main__":
    
    parents = dict()
    parents = {'N1': ['N2', 'N3', 'N4'], 'N3': ['N6', 'N7'], 'N4': ['N3'], 'N5': ['N4', 'N8'], 'N6': ['N13'],
               'N8': ['N9'], 'N9': ['N11'], 'N10': ['N7', 'N9'], 'N11': ['N14'], 'N12': ['N5']}

    print("Depth-first search:")
    dist = DFS_dist_from_node('N1', parents)
    print(dist)
    
    print("Breadth-first search:")
    dist = BFS_dist_from_node('N1', parents)
    print(dist)
    
