##################################
#                                #
#                                # 
#   Author: Usama Saqib          #
#   Description: UPGMA algorithm #
#                                #
##################################
import re
import copy

# enumerator for backtrack pointers
class Enum(set):
    def __getattr__(self, name):
        if name in self:
            return name
        raise AttributeError


# needleman wunsch implementation
def needleman_wunsch(seq1, seq2):
    match = 4
    mismatch = -3
    gap_penalty = -2

    BT = Enum(["LEFT", "UP", "DIAG_MATCH", "DIAG_MIS", "END"])

    seq1_len = len(seq1)
    seq2_len = len(seq2)


    # creates a matrix and initializes all cells to 0, and the backtrack pointers to End
    matrix = [ [ (0, BT.END) for x in range(seq2_len)] for y in range(seq1_len)]
    
    #calculate the score for each cell
    def score( ls1, ls2):

        #depending on the max value, set the backtrack pointer accordingly.
        def maxAndPointer(x, y, z, diag_type):
            m = max( x, y, z)
            if m == x:
                if (diag_type):
                    return (m, BT.DIAG_MATCH)
                else:
                    return (m, BT.DIAG_MIS)
            elif m == y:
                return (m, BT.LEFT)
            elif m == z:
                return (m, BT.UP)
        
        if ( seq1[ls1] == seq2[ls2]): 
            diag_match = True
            diag_score = matrix[ls1][ls2][0] + match
        else: 
            diag_match = False
            diag_score = matrix[ls1][ls2][0] + mismatch

        matrix[ls1][ls2] =  maxAndPointer( diag_score + matrix[ls1-1][ls2-1][0], (matrix[ls1][ls2-1])[0] + gap_penalty, (matrix[ls1-1][ls2])[0] + gap_penalty, diag_match)

    def total_distance( x, y, pointer, distance):
        if ( pointer == BT.END):
            return distance
        else:
            if (pointer == BT.LEFT):
                distance = total_distance( x-1, y, (matrix[y][x-1])[1], distance+1)
            elif ( pointer == BT.UP):  
                distance = total_distance( x, y-1, (matrix[y-1][x])[1], distance+1)
            elif (pointer == BT.DIAG_MIS):
                distance = total_distance( x-1, y-1, (matrix[y-1][x-1])[1], distance+1)
            elif (pointer == BT.DIAG_MATCH): 
                distance = total_distance( x-1, y-1, (matrix[y-1][x-1])[1], distance)

            return distance

    ### base conditions ###
    matrix[0] = [ (x * gap_penalty, BT.LEFT) if x > 0 else (0, BT.END) for x in range(seq2_len)]

    for y in range(1, seq1_len):
        matrix[y][0] = ( (matrix[y-1][0])[0] + gap_penalty, BT.UP) 

    ### matrix construction ###

    # runs score function for each cell value in the matrix.
    [[ score(y, x) for x in range(1, seq2_len)] for y in range(1, seq1_len)] 

    total_dist = total_distance( seq2_len-1, seq1_len-1, (matrix[seq1_len-1][seq2_len-1])[1], 0)

    return total_dist

def processFile():
    
    with open( "input.fa", 'r') as f:
        read_data = f.read()
    f.close()
    
    # returns all sequences in a list.
    expr = re.findall( r">.*\n([ATCG]+)", read_data)
    
    return expr

def upgma( distance_matrix, num_seqs):

    # the distance_matrix is repeadetly changed throughout the algorithm. A copy of the original
    # is kept in order to calculate distance between clusters.
    distance_matrix_original = copy.deepcopy(distance_matrix)

    # this list will be modified as the clusters are built. A cluster is represented by a tupule.
    # Example: (((0, 2), 3), (1, 4)) -> Let speciesX be denoted by sX.
    # This is a cluster of i:(s0 and s2), ii: (i) and s3, iii: s1 and s4, iv: (ii) and (iii)
    clusters = [ x for x in range( num_seqs) ]

    # this initializes a dictionary to store the distances between clusters.
    clusters_distances = { x : 0 for x in range(num_seqs) }    

    def findClosestCluster():
        smallest = 1000

        for key in distance_matrix:
            for val in distance_matrix[key]:
                if distance_matrix[key][val] < smallest:
                    smallest = distance_matrix[key][val]
                    cluster1 = key
                    cluster2 = val
         
        return (cluster1, cluster2)

    # implements the formula given during the lectures.
    def distFormula(c1, c2):
        
        # this return a generator.
        # the point of this function is to 'flatten' a nested tupule, so that
        # the necessary operations can be performed on it.
        # Example: ((0, (1, 3)), 2) -> [0, 1, 3, 2]
        def straighten(c):
            for i in c:
                if type(i) is tuple:
                    for j in straighten(i):
                        yield j
                else:
                    yield i

        numDist = 0
        lc1 = 0
        lc2 = 0

        flattened_c1 = straighten(c1)

        for leaf1 in flattened_c1:
            lc1 += 1
            lc2 = 0
            if ( type(c2) is tuple):
                for leaf2 in straighten(c2):
                    numDist += distance_matrix_original[min(leaf1, leaf2)][max(leaf1, leaf2)]
                    lc2 += 1
            else:
                numDist += distance_matrix_original[min(leaf1, c2)][max(leaf1, c2)]
                lc2 = 1

        denom = lc1 * lc2

        return numDist / denom

    # this function creates a new distance_matrix with the newly created cluster.
    def newDistanceMatrix( merge_clusters):

        # applies the distFormula to all other clusters.
        # this matrix is represented using a dictionary.
        new_dist = { x : distFormula(merge_clusters, x) for x in distance_matrix.keys() if x not in merge_clusters }

        del distance_matrix[merge_clusters[0]]
        del distance_matrix[merge_clusters[1]]

        # delete the individual clusters that are now merged.
        for key in distance_matrix:
            if merge_clusters[0] in distance_matrix[key].keys():
                del distance_matrix[key][merge_clusters[0]]
            if merge_clusters[1] in distance_matrix[key].keys():
                del distance_matrix[key][merge_clusters[1]]

        # create a new row in the matrix representing the new cluster
        distance_matrix[merge_clusters] = new_dist 

        return distance_matrix

    def output(c):
        out = ""
        if ( type(c) is int ):
            return "species " + str(c+1)
        else:
            out += "[" + output(c[0]) + ":" + str(clusters_distances[c[0]]) + "-" + output(c[1]) + ":" + str(clusters_distances[c[1]]) + "]"
            return out



    while len(clusters) > 1:
        merge_clusters = findClosestCluster()

        # remove the merged clusters
        clusters.remove(merge_clusters[0])
        clusters.remove(merge_clusters[1])

        new_height = distance_matrix[merge_clusters[0]][merge_clusters[1]] / 2

        distance_matrix = newDistanceMatrix(merge_clusters)

        clusters.append(merge_clusters)
 
        clusters_distances[merge_clusters[0]] = new_height - clusters_distances[merge_clusters[0]]
        clusters_distances[merge_clusters[1]] = new_height - clusters_distances[merge_clusters[1]]  
        clusters_distances[merge_clusters] = new_height #min(clusters_distances[merge_clusters[0]], clusters_distances[merge_clusters[1]]) 

    print(output(clusters[0]) + ":0.0")

def run():
    seqs = processFile()
    num_seqs = len(seqs)
    
    sequences = [ {"species-"+str(x) : "_"+seqs[x] } for x in range(num_seqs) ]

    # creates a new distance_matrix by calculating edit distance between all possible pairs of unique sequences, using the needleman-wunsch algorithm.
    distance_matrix = { y : { x : needleman_wunsch( sequences[y]["species-"+str(y)], sequences[x]["species-"+str(x)]) for x in range(y+1, num_seqs) } for y in range(num_seqs) }

    upgma( distance_matrix, num_seqs)


if __name__ == "__main__":
    run()
