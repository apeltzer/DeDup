package test.java;

import java.util.Set;
import java.util.stream.Stream;
import java.util.stream.Collectors;

public class Data {
       public static final Set<String> RMDupperTest__resolveDuplicate_yields_best_quality_expectedReadNames = Stream.of(
      "M_NS500559:30:HNM25BGXX:3:12601:6521:5315",
     	"M_NS500559:30:HNM25BGXX:1:12212:11266:3850",
     	"M_NS500559:30:HNM25BGXX:3:22507:7192:5163",
     	"M_NS500559:30:HNM25BGXX:4:21508:20245:5884",
     	"M_NS500559:30:HNM25BGXX:2:12204:22389:14660",
     	"M_NS500559:30:HNM25BGXX:2:11301:7589:8854",
     	"M_NS500559:30:HNM25BGXX:4:22605:13798:13816",
     	"M_NS500559:30:HNM25BGXX:1:13106:10220:17034",
     	"R_NS500559:30:HNM25BGXX:1:21203:15540:16052",
     	"M_NS500559:30:HNM25BGXX:1:13202:2193:2433",
     	"M_NS500559:30:HNM25BGXX:3:12607:11825:10419",
     	"M_NS500559:30:HNM25BGXX:4:23605:5727:2903",
     	"M_NS500559:30:HNM25BGXX:3:22512:6947:7877",
     	"M_NS500559:30:HNM25BGXX:1:11102:13689:7752",
     	"M_NS500559:30:HNM25BGXX:1:12111:13698:2649",
     	"M_NS500559:30:HNM25BGXX:1:12205:16021:11233",
     	"M_NS500559:30:HNM25BGXX:1:23307:13826:3081",
     	"M_NS500559:30:HNM25BGXX:3:21403:16465:9082",
     	"M_NS500559:30:HNM25BGXX:1:12304:7811:5239",
     	"M_NS500559:30:HNM25BGXX:1:23311:19533:3291",
     	"M_NS500559:30:HNM25BGXX:3:12504:16048:7514",
     	"M_NS500559:30:HNM25BGXX:4:22505:24153:2528").collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__resolveDuplicate_forward_expectedReadNames = Stream.of(
        "F_NS500559:30:HNM25BGXX:1:13209:14750:12764",
        "F_NS500559:30:HNM25BGXX:2:22309:10359:11690").collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__resolveDuplicate_forward_with_merged_expectedReadNames = Stream.of(
        "M_NS500559:30:HNM25BGXX:1:11102:13689:7752",
        "M_NS500559:30:HNM25BGXX:1:12111:13698:2649").collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__resolveDuplicate_reverse_expectedReadNames = Stream.of(
        "R_NS500559:30:HNM25BGXX:1:21203:15540:16052",
        "R_NS500559:30:HNM25BGXX:4:12601:1145:13392",
        "R_NS500559:30:HNM25BGXX:4:12406:10671:8628").collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__resolveDuplicate_reverse_with_merged_expectedReadNames = Stream.of(
        "M_NS500559:30:HNM25BGXX:1:12111:13698:2649",
        "M_NS500559:30:HNM25BGXX:1:11102:13689:7752").collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__resolveDuplicate_reverse_merged_expectedReadNames = Stream.of(
        "R_NS500559:30:HNM25BGXX:1:23312:9844:3571",
        "M_NS500559:30:HNM25BGXX:2:12204:22389:14660",
        "M_NS500559:30:HNM25BGXX:4:11506:2444:17383"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__queueOrOutput_nextAlignmentStartBeyondFirstEndTriggersCheckForDuplication_pre_expectedReadNames = Stream.of(
        "M_NS500559:30:HNM25BGXX:1:12108:16299:13403",
        "M_NS500559:30:HNM25BGXX:3:12601:6521:5315",
        "M_NS500559:30:HNM25BGXX:2:12204:22389:14659"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__queueOrOutput_nextAlignmentStartBeyondFirstEndTriggersCheckForDuplication_post_expectedReadNames = Stream.of(
        "M_NS500559:30:HNM25BGXX:3:12601:6521:5315",
        "M_NS500559:30:HNM25BGXX:2:12204:22389:14659",
        "M_NS500559:30:HNM25BGXX:2:12204:22389:14660"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__queueOrOutput_unsorted_read_one_forward__read_two_reverse_expectedReadNames = Stream.of(
        "F_1_1"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__queueOrOutput_unsorted_read_one_reverse__read_two_forward_expectedReadNames = Stream.of(
        "F_2_1"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__queueOrOutput_unsorted_read_two_forward__read_one_reverse_expectedReadNames = Stream.of(
        "R_3_1"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__queueOrOutput_unsorted_read_two_reverse__read_one_forward_expectedReadNames = Stream.of(
        "R_4_1"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest_allReadsAsMerged_expectedReadNames = Stream.of(
        "F_1",
        "M_1",
        "F_2"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__strandForward_expectedReadNames = Stream.of(
        "F_NS500559:30:HNM25BGXX:4:21508:20245:5884",
        "F_NS500559:30:HNM25BGXX:1:22107:23195:6187"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__strandReverse_expectedReadNames = Stream.of(
        "R_NS500559:30:HNM25BGXX:4:21508:20245:5884",
        "R_NS500559:30:HNM25BGXX:1:22107:23195:6187"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__single_expectedReadNames = Stream.of(
        "M_NS500559:30:HNM25BGXX:2:13301:8450:16718"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__stack_expectedReadNames = Stream.of(
        "M_DUPS_STACK",
        "M_DUPS_STACK_1_M_1",
        "M_DUPS_STACK_2_M_1",
        "M_DUPS_STACK_3_MFR_1",
        "M_UNIQ_STACK",
        "M_UNIQ_STACK_1_M_1",
        "M_UNIQ_STACK_1_M_3",
        "M_UNIQ_STACK_1_M_2",
        "M_UNIQ_STACK_2_M_1",
        "M_UNIQ_STACK_2_M_2",
        "M_UNIQ_STACK_2_M_3",
        "M_UNIQ_STACK_3_MFR_1",
        "F_UNIQ_STACK_3_MFR_1",
        "R_UNIQ_STACK_3_MFR_1"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__stackReadOne_expectedReadNames = Stream.of(
        "F_DUPS_STACK",
        "F_DUPS_STACK_1_M_1",
        "F_DUPS_STACK_2_M_1",
        "F_UNIQ_STACK",
        "F_UNIQ_STACK_1_M_1",
        "R_UNIQ_STACK_1_M_2",
        "F_UNIQ_STACK_1_M_3",
        "F_UNIQ_STACK_2_M_1",
        "F_UNIQ_STACK_2_M_2",
        "R_UNIQ_STACK_2_M_3"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__stackReadTwo_expectedReadNames = Stream.of(
        "R_DUPS_STACK",
        "R_DUPS_STACK_1_M_1",
        "R_DUPS_STACK_2_M_1",
        "R_UNIQ_STACK",
        "F_UNIQ_STACK_1_M_2",
        "R_UNIQ_STACK_1_M_1",
        "R_UNIQ_STACK_1_M_3",
        "R_UNIQ_STACK_2_M_1",
        "F_UNIQ_STACK_2_M_3",
        "R_UNIQ_STACK_2_M_2"
        ).collect(Collectors.toSet());
        public static final Set<String> RMDupperTest__MAPQZero_expectedReadNames = Stream.of(
        "F_NS500559:30:HNM25BGXX:4:21508:20245:5884",
        "F_mz_1"
        ).collect(Collectors.toSet());
}
