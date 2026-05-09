<?php
// ==================== ERROR DISPLAY ====================
ini_set('display_errors', 1);
ini_set('display_startup_errors', 1);
error_reporting(E_ALL);

header('Content-Type: text/html; charset=utf-8');
date_default_timezone_set('UTC');

// ==================== DB CONFIG ====================
$DB_HOST = 'xxx';
$DB_USER = 'xxx';
$DB_PASS = 'xxx';
$DB_NAME = 'xxx';

####  Provided table dumps as csv's in the directory for historical

// ==================== TABLES ====================
$LIQ_TABLE       = 'crossborder_capital_global_liquidity';
$MOVE_TABLE      = 'US-Treasury-Move-Index';
$BTC_TABLE       = 'Bitcoin-USD';
$RISKLOVE_TABLE  = 'Global-RiskLove-Composite';
$CHINA_CREDIT_IMPULSE_TABLE = 'China-Bloomberg-Credit-Impulse-Index';

// ==================== CHINA CREDIT IMPULSE SETTINGS ====================
// Change this one value to control the China Credit Impulse impact lead/lag everywhere relevant:
// - blended wave calculation
// - standalone China Credit Impulse plotted trace
// - standalone China Credit Impulse deviation bands
// Positive value shifts the impact/plot forward by N months.
$CHINA_CREDIT_IMPULSE_SHIFT_MONTHS = 9;
$CHINA_CREDIT_IMPULSE_DOMINANCE_MULTIPLE = 0.65;

// ==================== CONNECT ====================
$mysqli = new mysqli($DB_HOST, $DB_USER, $DB_PASS, $DB_NAME);
if ($mysqli->connect_errno) {
    die("DB connection failed: " . $mysqli->connect_error);
}

// ==================== HELPERS ====================
function dmy_to_ymd($dmy){
    $dt = DateTime::createFromFormat('d/m/Y',$dmy);
    return $dt ? $dt->format('Y-m-d') : null;
}
function mean($a){
    $a = array_values(array_filter($a, fn($v) => $v !== null && is_numeric($v)));
    return count($a) ? array_sum($a)/count($a) : 0;
}
function stddev($a){
    $a = array_values(array_filter($a, fn($v) => $v !== null && is_numeric($v)));
    if(count($a) < 1) return 0;
    $m=mean($a); $v=0;
    foreach($a as $x){ $v+=pow($x-$m,2); }
    return sqrt($v/count($a));
}
function rolling_median($a,$w){
    $o=[];
    for($i=0;$i<count($a);$i++){
        if($i<$w){ $o[]=null; continue; }
        $s=array_slice($a,$i-$w,$w);
        $s=array_values(array_filter($s, fn($v)=>$v!==null && is_numeric($v)));
        if(!count($s)){ $o[]=null; continue; }
        sort($s);
        $o[]=$s[floor(count($s)/2)];
    }
    return $o;
}
function ema($v,$s){
    $o=[]; $k=2/($s+1); $e=null;
    foreach($v as $x){
        if($x===null){
            $o[]=$e;
            continue;
        }
        $e=($e===null)?$x:(($x-$e)*$k+$e);
        $o[]=$e;
    }
    return $o;
}
function normalize($a){
    $clean = array_values(array_filter($a, fn($v)=>$v!==null && is_numeric($v)));
    $m = mean($clean);
    $s = stddev($clean);
    return array_map(fn($x)=>($x!==null && $s>0?($x-$m)/$s:null), $a);
}
function shift_array($arr, $shift){
    $out = array_fill(0, count($arr), null);
    for($i=0;$i<count($arr);$i++){
        $j = $i + $shift;
        if($j >= 0 && $j < count($arr)) $out[$j] = $arr[$i];
    }
    return $out;
}
function add_months_ymd($date, $months){
    $dt  = new DateTime($date);
    $day = (int)$dt->format('d');
    $dt->modify(($months >= 0 ? '+' : '') . $months . ' months');
    if ((int)$dt->format('d') !== $day) {
        $dt->modify('last day of previous month');
    }
    return $dt->format('Y-m-d');
}
function build_shifted_map($dates, $vals, $shiftMonths, $fromDate = null){
    $out = [];
    $n = min(count($dates), count($vals));
    for($i=0; $i<$n; $i++){
        if($vals[$i] === null) continue;
        $d = add_months_ymd($dates[$i], $shiftMonths);
        if($fromDate !== null && $d < $fromDate) continue;
        $out[$d] = $vals[$i];
    }
    ksort($out);
    return $out;
}
function align_two_series($a, $b, $fromDate = null){
    $dates = array_values(array_intersect(array_keys($a), array_keys($b)));
    sort($dates);

    $x = []; $y = []; $d = [];
    foreach($dates as $dt){
        if($fromDate !== null && $dt < $fromDate) continue;
        if($a[$dt] === null || $b[$dt] === null) continue;
        $d[] = $dt;
        $x[] = (float)$a[$dt];
        $y[] = (float)$b[$dt];
    }
    return [$d, $x, $y];
}
function pearson_corr($x, $y){
    $pairs = [];
    $n = min(count($x), count($y));
    for($i=0; $i<$n; $i++){
        if($x[$i]===null || $y[$i]===null) continue;
        $pairs[] = [(float)$x[$i], (float)$y[$i]];
    }
    if(count($pairs) < 3) return null;

    $xs = array_column($pairs, 0);
    $ys = array_column($pairs, 1);
    $mx = mean($xs);
    $my = mean($ys);
    $num = 0.0; $dx = 0.0; $dy = 0.0;

    for($i=0; $i<count($pairs); $i++){
        $a = $xs[$i] - $mx;
        $b = $ys[$i] - $my;
        $num += $a * $b;
        $dx  += $a * $a;
        $dy  += $b * $b;
    }

    if($dx <= 0 || $dy <= 0) return null;
    return $num / sqrt($dx * $dy);
}
function r_squared_from_corr($corr){
    return ($corr === null) ? null : $corr * $corr;
}
function directional_hit_rate($x, $y){
    $n = min(count($x), count($y));
    if($n < 4) return null;

    $hits = 0; $tot = 0;
    for($i=1; $i<$n; $i++){
        if($x[$i]===null || $x[$i-1]===null || $y[$i]===null || $y[$i-1]===null) continue;
        $dx = $x[$i] - $x[$i-1];
        $dy = $y[$i] - $y[$i-1];
        if($dx == 0 || $dy == 0) continue;
        $hits += ((($dx > 0) && ($dy > 0)) || (($dx < 0) && ($dy < 0))) ? 1 : 0;
        $tot++;
    }
    return $tot ? $hits / $tot : null;
}
function turning_point_hit_rate($x, $y){
    $n = min(count($x), count($y));
    if($n < 6) return null;

    $hits = 0; $tot = 0;
    for($i=2; $i<$n; $i++){
        if(
            $x[$i]===null || $x[$i-1]===null || $x[$i-2]===null ||
            $y[$i]===null || $y[$i-1]===null || $y[$i-2]===null
        ) continue;

        $dx1 = $x[$i-1] - $x[$i-2];
        $dx2 = $x[$i]   - $x[$i-1];
        $dy1 = $y[$i-1] - $y[$i-2];
        $dy2 = $y[$i]   - $y[$i-1];

        $turnX = ($dx1 > 0 && $dx2 < 0) || ($dx1 < 0 && $dx2 > 0);
        $turnY = ($dy1 > 0 && $dy2 < 0) || ($dy1 < 0 && $dy2 > 0);

        if(!$turnX) continue;
        $hits += $turnY ? 1 : 0;
        $tot++;
    }
    return $tot ? $hits / $tot : null;
}
function composite_impact_score($corr, $r2, $dirHit, $turnHit){
    $corr    = ($corr    === null) ? 0 : abs($corr);
    $r2      = ($r2      === null) ? 0 : $r2;
    $dirHit  = ($dirHit  === null) ? 0 : $dirHit;
    $turnHit = ($turnHit === null) ? 0 : $turnHit;

    return 100 * (
        0.40 * $corr +
        0.25 * $r2 +
        0.20 * $dirHit +
        0.15 * $turnHit
    );
}
function pct_rank_latest($arr){
    $clean = array_values(array_filter($arr, fn($v) => $v !== null && is_numeric($v)));
    if(!$clean) return null;
    $latest = end($clean);
    sort($clean);
    $count = count($clean);
    $le = 0;
    foreach($clean as $v){
        if($v <= $latest) $le++;
    }
    return $count ? $le / $count : null;
}
function simple_ma($arr, $window){
    $out = [];
    for($i=0; $i<count($arr); $i++){
        if($i < $window - 1){
            $out[] = null;
            continue;
        }
        $slice = array_slice($arr, $i - $window + 1, $window);
        $slice = array_values(array_filter($slice, fn($v)=>$v!==null));
        $out[] = count($slice) ? mean($slice) : null;
    }
    return $out;
}
function months_between($d1, $d2){
    $a = new DateTime($d1);
    $b = new DateTime($d2);
    $y = ((int)$b->format('Y') - (int)$a->format('Y')) * 12;
    $m = (int)$b->format('n') - (int)$a->format('n');
    return $y + $m;
}
function safe_max_non_null($arr){
    $arr = array_values(array_filter($arr, fn($v)=>$v!==null && is_numeric($v)));
    return count($arr) ? max($arr) : null;
}
function safe_min_non_null($arr){
    $arr = array_values(array_filter($arr, fn($v)=>$v!==null && is_numeric($v)));
    return count($arr) ? min($arr) : null;
}
function signed_log_value($v, $scale = 1.0){
    if($v === null || !is_numeric($v)) return null;
    $v = (float)$v;
    if($v == 0.0) return 0.0;
    return ($v < 0 ? -1 : 1) * log(1 + abs($v) * $scale);
}
function signed_log_series($arr, $scale = 1.0){
    return array_map(fn($v) => signed_log_value($v, $scale), $arr);
}
function series_to_map($dates, $vals, $fromDate = null){
    $out = [];
    $n = min(count($dates), count($vals));
    for($i=0; $i<$n; $i++) {
        if($dates[$i] === null || $vals[$i] === null) continue;
        if($fromDate !== null && $dates[$i] < $fromDate) continue;
        $out[$dates[$i]] = (float)$vals[$i];
    }
    ksort($out);
    return $out;
}
function compress_regimes($dates, $states){
    $out = [];
    if(!$dates || !$states || count($dates) !== count($states)) return $out;

    $start = $dates[0];
    $state = $states[0];

    for($i=1; $i<count($dates); $i++){
        if($states[$i] !== $state){
            $out[] = ['start'=>$start, 'end'=>$dates[$i-1], 'state'=>$state];
            $start = $dates[$i];
            $state = $states[$i];
        }
    }
    $out[] = ['start'=>$start, 'end'=>$dates[count($dates)-1], 'state'=>$state];
    return $out;
}

// ==================== PEAK / BOTTOM FINDERS ====================
function find_btc_cycle_peaks($btcMonthlyMap, $lookback = 5, $lookforward = 5, $futureWindow = 12, $minDrawdown = 0.25, $minSepMonths = 8){
    $keys = array_keys($btcMonthlyMap);
    sort($keys);

    $peaks = [];
    $lastPeak = null;

    for($i=$lookback; $i<count($keys)-$lookforward; $i++){
        $d = $keys[$i];
        $v = $btcMonthlyMap[$d];

        $prev = [];
        for($j=$i-$lookback; $j<$i; $j++) $prev[] = $btcMonthlyMap[$keys[$j]];
        $next = [];
        for($j=$i+1; $j<=$i+$lookforward; $j++) $next[] = $btcMonthlyMap[$keys[$j]];

        if($v <= max($prev) || $v < max($next)) continue;

        $futureMin = $v;
        $end = min(count($keys)-1, $i+$futureWindow);
        for($j=$i+1; $j<=$end; $j++){
            $futureMin = min($futureMin, $btcMonthlyMap[$keys[$j]]);
        }

        $dd = ($futureMin / $v) - 1.0;
        if($dd > -$minDrawdown) continue;

        if($lastPeak !== null && months_between($lastPeak, $d) < $minSepMonths) continue;

        $peaks[] = $d;
        $lastPeak = $d;
    }

    return $peaks;
}
function find_btc_cycle_bottoms($btcMonthlyMap, $lookback = 5, $lookforward = 5, $futureWindow = 12, $minRally = 0.35, $minSepMonths = 8){
    $keys = array_keys($btcMonthlyMap);
    sort($keys);

    $bottoms = [];
    $lastBottom = null;

    for($i=$lookback; $i<count($keys)-$lookforward; $i++){
        $d = $keys[$i];
        $v = $btcMonthlyMap[$d];

        $prev = [];
        for($j=$i-$lookback; $j<$i; $j++) $prev[] = $btcMonthlyMap[$keys[$j]];
        $next = [];
        for($j=$i+1; $j<=$i+$lookforward; $j++) $next[] = $btcMonthlyMap[$keys[$j]];

        if($v >= min($prev) || $v > min($next)) continue;

        $futureMax = $v;
        $end = min(count($keys)-1, $i+$futureWindow);
        for($j=$i+1; $j<=$end; $j++){
            $futureMax = max($futureMax, $btcMonthlyMap[$keys[$j]]);
        }

        $rally = ($futureMax / $v) - 1.0;
        if($rally < $minRally) continue;

        if($lastBottom !== null && months_between($lastBottom, $d) < $minSepMonths) continue;

        $bottoms[] = $d;
        $lastBottom = $d;
    }

    return $bottoms;
}

// ==================== CONDITION PACKS ====================
function condition_pack_peak($dates, $btcVals, $waveVals){
    $n = min(count($dates), count($btcVals), count($waveVals));
    if($n < 7){
        return [
            'btc_hot' => false,
            'wave_rollover_3m' => false,
            'wave_below_3m_high' => false,
            'wave_below_6m_ma' => false,
            'btc_up_wave_down_2m' => false,
            'corr6_lt_corr12' => false
        ];
    }

    $btcPct = pct_rank_latest($btcVals);
    $ma6 = simple_ma($waveVals, 6);

    $corr6  = pearson_corr(array_slice($btcVals, -6),  array_slice($waveVals, -6));
    $corr12 = pearson_corr(array_slice($btcVals, -12), array_slice($waveVals, -12));

    $latestWave = $waveVals[$n-1] ?? null;
    $wave3Ago   = $waveVals[$n-4] ?? null;
    $wave2Ago   = $waveVals[$n-3] ?? null;
    $wave3mHigh = safe_max_non_null(array_slice($waveVals, -3));
    $waveMa6    = $ma6[$n-1] ?? null;

    $latestBTC  = $btcVals[$n-1] ?? null;
    $btc2Ago    = $btcVals[$n-3] ?? null;

    return [
        'btc_hot'              => ($btcPct !== null && $btcPct >= 0.75),
        'wave_rollover_3m'     => ($wave3Ago !== null && $latestWave !== null && $latestWave < $wave3Ago),
        'wave_below_3m_high'   => ($latestWave !== null && $wave3mHigh !== null && $latestWave < $wave3mHigh),
        'wave_below_6m_ma'     => ($waveMa6 !== null && $latestWave !== null && $latestWave < $waveMa6),
        'btc_up_wave_down_2m'  => ($btc2Ago !== null && $wave2Ago !== null && $latestBTC !== null && $latestWave !== null && ($latestBTC > $btc2Ago) && ($latestWave < $wave2Ago)),
        'corr6_lt_corr12'      => ($corr6 !== null && $corr12 !== null && $corr6 < $corr12)
    ];
}
function condition_pack_bottom($dates, $btcVals, $waveVals){
    $n = min(count($dates), count($btcVals), count($waveVals));
    if($n < 7){
        return [
            'btc_cold' => false,
            'wave_turnup_3m' => false,
            'wave_above_3m_low' => false,
            'wave_above_6m_ma' => false,
            'btc_down_wave_up_2m' => false,
            'corr6_gt_corr12' => false
        ];
    }

    $btcPct = pct_rank_latest($btcVals);
    $ma6 = simple_ma($waveVals, 6);

    $corr6  = pearson_corr(array_slice($btcVals, -6),  array_slice($waveVals, -6));
    $corr12 = pearson_corr(array_slice($btcVals, -12), array_slice($waveVals, -12));

    $latestWave = $waveVals[$n-1] ?? null;
    $wave3Ago   = $waveVals[$n-4] ?? null;
    $wave2Ago   = $waveVals[$n-3] ?? null;
    $wave3mLow  = safe_min_non_null(array_slice($waveVals, -3));
    $waveMa6    = $ma6[$n-1] ?? null;

    $latestBTC  = $btcVals[$n-1] ?? null;
    $btc2Ago    = $btcVals[$n-3] ?? null;

    return [
        'btc_cold'             => ($btcPct !== null && $btcPct <= 0.25),
        'wave_turnup_3m'       => ($wave3Ago !== null && $latestWave !== null && $latestWave > $wave3Ago),
        'wave_above_3m_low'    => ($latestWave !== null && $wave3mLow !== null && $latestWave > $wave3mLow),
        'wave_above_6m_ma'     => ($waveMa6 !== null && $latestWave !== null && $latestWave > $waveMa6),
        'btc_down_wave_up_2m'  => ($btc2Ago !== null && $wave2Ago !== null && $latestBTC !== null && $latestWave !== null && ($latestBTC < $btc2Ago) && ($latestWave > $wave2Ago)),
        'corr6_gt_corr12'      => ($corr6 !== null && $corr12 !== null && $corr6 > $corr12)
    ];
}

// ==================== DRIVER SCORERS ====================
function score_peak_driver_for_wave($peakDatePlot, $fullDates, $btcSeries, $waveSeries){
    $idx = array_search($peakDatePlot, $fullDates, true);
    if($idx === false || $idx < 12) return null;

    $dates = array_slice($fullDates, max(0, $idx-11), 12);
    $btc   = array_slice($btcSeries, max(0, $idx-11), 12);
    $wave  = array_slice($waveSeries, max(0, $idx-11), 12);

    $corr12 = pearson_corr($btc, $wave);

    $wavePeakDate = null;
    $wavePeakVal  = null;
    for($i=0; $i<count($wave); $i++){
        if($wave[$i] === null) continue;
        if($wavePeakVal === null || $wave[$i] > $wavePeakVal){
            $wavePeakVal  = $wave[$i];
            $wavePeakDate = $dates[$i];
        }
    }
    if($wavePeakDate === null){
        return null;
    }

    $leadMonths = months_between($wavePeakDate, $peakDatePlot);
    $leadScore  = 0.0;
    if($leadMonths >= 0 && $leadMonths <= 6){
        $leadScore = 1 - abs($leadMonths - 2) / 6;
        if($leadScore < 0) $leadScore = 0;
    }

    $rollover = false;
    if(count($wave) >= 4 && $wave[count($wave)-1] !== null && $wave[count($wave)-4] !== null){
        $rollover = ($wave[count($wave)-1] < $wave[count($wave)-4]);
    }

    $peakScore = 100 * (
        0.55 * max(0, (float)$corr12) +
        0.30 * $leadScore +
        0.15 * ($rollover ? 1 : 0)
    );

    return [
        'corr12'     => $corr12,
        'leadMonths' => $leadMonths,
        'rollover'   => $rollover,
        'peakScore'  => $peakScore
    ];
}
function score_bottom_driver_for_wave($bottomDatePlot, $fullDates, $btcSeries, $waveSeries){
    $idx = array_search($bottomDatePlot, $fullDates, true);
    if($idx === false || $idx < 12) return null;

    $dates = array_slice($fullDates, max(0, $idx-11), 12);
    $btc   = array_slice($btcSeries, max(0, $idx-11), 12);
    $wave  = array_slice($waveSeries, max(0, $idx-11), 12);

    $corr12 = pearson_corr($btc, $wave);

    $waveBottomDate = null;
    $waveBottomVal  = null;
    for($i=0; $i<count($wave); $i++){
        if($wave[$i] === null) continue;
        if($waveBottomVal === null || $wave[$i] < $waveBottomVal){
            $waveBottomVal  = $wave[$i];
            $waveBottomDate = $dates[$i];
        }
    }
    if($waveBottomDate === null){
        return null;
    }

    $leadMonths = months_between($waveBottomDate, $bottomDatePlot);
    $leadScore  = 0.0;
    if($leadMonths >= 0 && $leadMonths <= 6){
        $leadScore = 1 - abs($leadMonths - 2) / 6;
        if($leadScore < 0) $leadScore = 0;
    }

    $turnup = false;
    if(count($wave) >= 4 && $wave[count($wave)-1] !== null && $wave[count($wave)-4] !== null){
        $turnup = ($wave[count($wave)-1] > $wave[count($wave)-4]);
    }

    $bottomScore = 100 * (
        0.55 * max(0, (float)$corr12) +
        0.30 * $leadScore +
        0.15 * ($turnup ? 1 : 0)
    );

    return [
        'corr12'      => $corr12,
        'leadMonths'  => $leadMonths,
        'turnup'      => $turnup,
        'bottomScore' => $bottomScore
    ];
}
function evaluate_template_match($current, $rates, $minRate = 0.60){
    $need = 0; $met = 0; $active = [];
    foreach($rates as $k => $rate){
        if($rate >= $minRate){
            $need++;
            $isMet = !empty($current[$k]);
            if($isMet) $met++;
            $active[$k] = [
                'historicalRate' => $rate,
                'currentMet' => $isMet
            ];
        }
    }

    $status = 'No';
    if($need > 0){
        $ratio = $met / $need;
        if($ratio >= 0.80)      $status = 'Yes';
        elseif($ratio >= 0.50)  $status = 'Partial';
        else                    $status = 'No';
    }

    return [
        'need'   => $need,
        'met'    => $met,
        'status' => $status,
        'items'  => $active
    ];
}

// ==================== BLENDED WAVE WEIGHTS ====================
$W_LIQ = 0.50;
$W_QE  = 0.30;
$W_YCC = 0.20;

// ==================== PLOT START (GLOBAL) ====================
$PLOT_START = '2019-08-01';

// ==================== LOAD LIQUIDITY (PRIMARY) ====================
$liq_stock=[];
$q="
SELECT date_period, global_liquidity_tn$, shadow_monetary_base_tn$
FROM `$LIQ_TABLE`
WHERE global_liquidity_tn$ IS NOT NULL
  AND shadow_monetary_base_tn$ IS NOT NULL
ORDER BY STR_TO_DATE(date_period,'%d/%m/%Y') ASC";
$r=$mysqli->query($q);
while($row=$r->fetch_assoc()){
    $d=dmy_to_ymd($row['date_period']);
    if(!$d) continue;
    $liq_stock[$d]=floatval($row['global_liquidity_tn$'])
                  +floatval($row['shadow_monetary_base_tn$']);
}
$liq_dates=array_keys($liq_stock);
$liq_vals =array_values($liq_stock);

// ==================== YOY LOG IMPULSE ====================
$imp_dates=[]; $imp_vals=[];
for($i=12;$i<count($liq_vals);$i++){
    if($liq_vals[$i-12] <= 0 || $liq_vals[$i] <= 0) continue;
    $imp_dates[]=$liq_dates[$i];
    $imp_vals[]=log($liq_vals[$i]/$liq_vals[$i-12]);
}

// ==================== LOAD MOVE (LEVEL) ====================
$move_raw=[];
$r=$mysqli->query("SELECT date,value FROM `$MOVE_TABLE` ORDER BY date ASC");
while($row=$r->fetch_assoc()){
    $move_raw[$row['date']]=floatval($row['value']);
}
$move_vals=[];
foreach($imp_dates as $d){ $move_vals[]=$move_raw[$d]??null; }

// ==================== MOVE SUPPRESSION ====================
$move_smooth=rolling_median($move_vals,26);
$mc=array_filter($move_smooth,fn($v)=>$v!==null);
$mm=mean($mc); $ms=stddev($mc);

$vol_damp=[];
foreach($move_smooth as $v){
    if($v===null){ $vol_damp[]=null; continue; }
    $z=($ms > 0) ? (($v-$mm)/$ms) : 0;
    $d=exp(-0.25*$z);
    $vol_damp[]=max(0.6,min(1.4,$d));
}

// ==================== APPLY DAMPING ====================
$imp_adj=[]; $adj_dates=[];
for($i=0;$i<count($imp_vals);$i++){
    if($vol_damp[$i]===null) continue;
    $imp_adj[]=$imp_vals[$i]*$vol_damp[$i];
    $adj_dates[]=$imp_dates[$i];
}

// ==================== PRIMARY LIQUIDITY WAVE ====================
$w1=ema($imp_adj,4);
$w2=ema($w1,6);
$wave_raw=ema($w2,4);
$wave_z = normalize($wave_raw);

$acc=[];
for($i=1;$i<count($wave_z);$i++) {
    $acc[] = ($wave_z[$i]!==null && $wave_z[$i-1]!==null) ? ($wave_z[$i]-$wave_z[$i-1]) : null;
}
$acc_s=ema($acc,1);

$liq_wave_final=[isset($wave_z[0]) ? $wave_z[0] : null];
for($i=1;$i<count($wave_z);$i++){
    $liq_wave_final[] = ($wave_z[$i]!==null && isset($acc_s[$i-1]) && $acc_s[$i-1]!==null)
        ? ($wave_z[$i]+0.9*$acc_s[$i-1])
        : $wave_z[$i];
}

// ==================== LOAD YCC & QE ====================
function fetch_txt($u){
    $l=@file($u,FILE_IGNORE_NEW_LINES|FILE_SKIP_EMPTY_LINES);
    $o=[];
    if(!$l) return $o;
    foreach($l as $x){
        $p=explode(',',$x);
        if(count($p)==2 && is_numeric($p[1])) $o[$p[0]]=floatval($p[1]);
    }
    return $o;
}
function carry_forward($dates,$src){
    $o=[];$last=null;
    if(!is_array($src) || !count($src)){
        foreach($dates as $_) $o[]=null;
        return $o;
    }
    ksort($src);
    $keys=array_keys($src);
    $idx=0;
    $n=count($keys);
    foreach($dates as $d){
        while($idx < $n && $keys[$idx] <= $d){
            $last=$src[$keys[$idx]];
            $idx++;
        }
        $o[]=$last;
    }
    return $o;
}

// ==================== LOAD MONTHLY DB SERIES FOR STEALTH-QE STYLE WAVES ====================
function load_monthly_db_series_map($mysqli, $table){
    $out=[];
    $tableEsc=$mysqli->real_escape_string($table);
    $sql="SELECT date,value FROM `".$tableEsc."` WHERE date IS NOT NULL AND value IS NOT NULL ORDER BY date ASC";
    $res=$mysqli->query($sql);
    if($res){
        while($row=$res->fetch_assoc()){
            if(!isset($row['date'], $row['value']) || !is_numeric($row['value'])) continue;
            $ts=strtotime((string)$row['date']);
            if($ts===false) continue;
            $d=date('Y-m-01', $ts);
            $out[$d]=(float)$row['value'];
        }
    }
    ksort($out);
    return $out;
}


// ==================== DEBT-LIQUIDITY FORWARD ENGINE ====================
// Adds rollover + interest + Fed/Treasury liquidity requirement projection into the existing wave stack.
$FRED_API_KEY = 'xxxx'; ##USE FROM FRED (FREE)
$DEBT_LIQ_FORWARD_YEARS = 5;
$DEBT_LIQ_REQUIRED_MULTIPLE = 2.70;
$DEBT_LIQ_ROLLOVER_CYCLES = 8;
$DEBT_LIQ_CACHE_TTL = 86400;
$DEBT_LIQ_CACHE_DIR = __DIR__ . '/cache';
if(!is_dir($DEBT_LIQ_CACHE_DIR)) @mkdir($DEBT_LIQ_CACHE_DIR, 0755, true);

function dlq_cache_get($file, $ttl){
    if(file_exists($file) && (time() - filemtime($file)) < $ttl){
        $j = @file_get_contents($file);
        $d = json_decode($j, true);
        if(is_array($d)) return $d;
    }
    return null;
}
function dlq_cache_set($file, $data){ @file_put_contents($file, json_encode($data)); }
function dlq_fetch_json($url, $timeout = 30){
    $ctx = stream_context_create(['http'=>['method'=>'GET','header'=>"User-Agent: Mozilla/5.0\r\n",'timeout'=>$timeout]]);
    $json = @file_get_contents($url, false, $ctx);
    if($json === false || trim($json)==='') return null;
    $decoded = json_decode($json, true);
    return is_array($decoded) ? $decoded : null;
}
function dlq_fred_series_millions_to_tn($seriesId, $apiKey, $startDate, $endDate, $cacheDir, $ttl){
    $cache = $cacheDir . '/fred_' . strtolower($seriesId) . '_debt_liq_cache.json';
    $data = dlq_cache_get($cache, $ttl);
    if($data === null){
        $url = 'https://api.stlouisfed.org/fred/series/observations?' . http_build_query([
            'series_id'=>$seriesId,
            'api_key'=>$apiKey,
            'file_type'=>'json',
            'observation_start'=>$startDate,
            'observation_end'=>$endDate
        ]);
        $data = dlq_fetch_json($url, 30);
        if($data && isset($data['observations'])) dlq_cache_set($cache, $data);
    }
    $out = [];
    if($data && isset($data['observations']) && is_array($data['observations'])){
        foreach($data['observations'] as $obs){
            if(!isset($obs['date'], $obs['value']) || $obs['value']==='.' || !is_numeric($obs['value'])) continue;
            $out[$obs['date']] = ((float)$obs['value']) / 1000000.0; // FRED millions -> trillions
        }
    }
    ksort($out);
    return $out;
}
function dlq_month_start($date){ return date('Y-m-01', strtotime($date)); }
function dlq_monthly_dates($start, $end){
    $out=[];
    $dt = new DateTime(dlq_month_start($start));
    $stop = new DateTime(dlq_month_start($end));
    while($dt <= $stop){ $out[] = $dt->format('Y-m-d'); $dt->modify('+1 month'); }
    return $out;
}
function dlq_carry_map($dates, $src){
    $out=[]; $last=null; ksort($src);
    $srcKeys=array_keys($src); $j=0; $n=count($srcKeys);
    foreach($dates as $d){
        while($j<$n && $srcKeys[$j] <= $d){ $last=$src[$srcKeys[$j]]; $j++; }
        $out[$d]=$last;
    }
    return $out;
}
function dlq_normalize_amount_to_b($v){
    $v = str_replace([',','$'], '', trim((string)$v));
    return is_numeric($v) ? ((float)$v / 1000000000.0) : 0.0;
}
function dlq_classify_cluster($type, $term, $cmb){
    $type=strtoupper(trim((string)$type)); $term=strtoupper(trim((string)$term)); $cmb=strtoupper(trim((string)$cmb));
    if($type==='BILL' || $cmb==='YES') return 'Bills <1Y';
    if(strpos($term,'2-YEAR')!==false || strpos($term,'3-YEAR')!==false) return '2-3Y';
    if(strpos($term,'5-YEAR')!==false || strpos($term,'7-YEAR')!==false) return '5-7Y';
    if(strpos($term,'10-YEAR')!==false) return '10Y';
    if(strpos($term,'20-YEAR')!==false || strpos($term,'30-YEAR')!==false) return '20-30Y';
    if($type==='FRN') return 'FRN';
    if($type==='TIPS' || strpos($term,'TIPS')!==false) return 'TIPS';
    return 'Other';
}
function dlq_fetch_db_sum($mysqli, $table, $startDate){
    $sum=0.0;
    $sql="SELECT date,value FROM `".$mysqli->real_escape_string($table)."` WHERE date >= ? AND value IS NOT NULL ORDER BY date ASC";
    $stmt=$mysqli->prepare($sql);
    if(!$stmt) return 0.0;
    $stmt->bind_param('s', $startDate);
    $stmt->execute();
    $res=$stmt->get_result();
    while($row=$res->fetch_assoc()) if(isset($row['value']) && is_numeric($row['value'])) $sum += (float)$row['value'];
    $stmt->close();
    return $sum;
}
function dlq_dynamic_rollover_weights($bills, $notes, $bonds){
    $total=max(0,$bills)+max(0,$notes)+max(0,$bonds);
    if($total<=0){
        return ['Bills <1Y'=>0.40,'2-3Y'=>0.20,'5-7Y'=>0.25,'10Y'=>0.10,'20-30Y'=>0.05,'FRN'=>0,'TIPS'=>0,'Other'=>0];
    }
    return ['Bills <1Y'=>$bills/$total,'2-3Y'=>($notes/$total)*0.40,'5-7Y'=>($notes/$total)*0.40,'10Y'=>($notes/$total)*0.20,'20-30Y'=>$bonds/$total,'FRN'=>0,'TIPS'=>0,'Other'=>0];
}
function dlq_fetch_treasury_maturities_by_year($startYear, $endYear, $cacheDir, $ttl){
    $cache=$cacheDir.'/treasury_debt_liq_maturities_cache.json';
    $cached=dlq_cache_get($cache, $ttl);
    if(is_array($cached)) return $cached;
    $api='https://api.fiscaldata.treasury.gov/services/api/fiscal_service/v1/accounting/od/auctions_query';
    $page=1; $pageSize=10000; $years=[];
    do{
        $url=$api.'?'.http_build_query([
            'fields'=>'record_date,auction_date,issue_date,maturity_date,security_type,security_term,total_accepted,offering_amt,cash_management_bill_cmb',
            'filter'=>'maturity_date:gte:'.$startYear.'-01-01',
            'page[size]'=>$pageSize,
            'page[number]'=>$page,
            'sort'=>'maturity_date'
        ]);
        $d=dlq_fetch_json($url, 30);
        if(!$d || !isset($d['data']) || !is_array($d['data'])) break;
        foreach($d['data'] as $row){
            if(empty($row['maturity_date'])) continue;
            $y=(int)date('Y', strtotime($row['maturity_date']));
            if($y<$startYear || $y>$endYear) continue;
            $amount=0.0;
            if(isset($row['total_accepted'])) $amount=dlq_normalize_amount_to_b($row['total_accepted']);
            if($amount<=0 && isset($row['offering_amt'])) $amount=dlq_normalize_amount_to_b($row['offering_amt']);
            if($amount<=0) continue;
            $cluster=dlq_classify_cluster($row['security_type'] ?? '', $row['security_term'] ?? '', $row['cash_management_bill_cmb'] ?? '');
            if(!isset($years[$y])) $years[$y]=[];
            if(!isset($years[$y][$cluster])) $years[$y][$cluster]=0.0;
            $years[$y][$cluster]+=$amount;
        }
        $count=count($d['data']); $page++;
    }while($count===$pageSize && $page<100);
    dlq_cache_set($cache, $years);
    return $years;
}
function dlq_project_rollovers_annual_tn($actualClusters, $currentYear, $endYear, $weights, $cycles){
    $terms=['Bills <1Y'=>1,'2-3Y'=>3,'5-7Y'=>5,'10Y'=>10,'20-30Y'=>20,'FRN'=>2,'TIPS'=>10,'Other'=>5];
    $events=[]; $out=[];
    foreach($actualClusters as $year=>$clusters){
        if((int)$year < $currentYear) continue;
        foreach($clusters as $cluster=>$amountB){ if($amountB>0) $events[]=['year'=>(int)$year,'amountB'=>(float)$amountB]; }
    }
    for($cycle=1;$cycle<=$cycles;$cycle++){
        $next=[];
        foreach($events as $ev){
            foreach($weights as $dest=>$w){
                if($w<=0) continue;
                $fy=$ev['year'] + ($terms[$dest] ?? 5);
                if($fy>$endYear) continue;
                $amtB=$ev['amountB']*$w;
                if(!isset($out[$fy])) $out[$fy]=0.0;
                $out[$fy]+=$amtB/1000.0; // B -> Tn
                $next[]=['year'=>$fy,'amountB'=>$amtB];
            }
        }
        $events=$next;
        if(empty($events)) break;
    }
    ksort($out);
    return $out;
}
function dlq_latest_value($arr){ $v=null; foreach($arr as $x){ if($x!==null && is_numeric($x)) $v=$x; } return $v; }
function dlq_last_date_in_map($arr){ $last=null; foreach($arr as $d=>$v){ if($v!==null && is_numeric($v)) $last=$d; } return $last; }
function dlq_map_no_carry($dates, $src){ $out=[]; foreach($dates as $d){ $out[] = array_key_exists($d, $src) ? $src[$d] : null; } return $out; }

$debt_liq_future_end = date('Y-m-d', strtotime('+' . $DEBT_LIQ_FORWARD_YEARS . ' years'));

// Extend the core monthly timeline dynamically so the wave traces and sigma/deviation bands
// continue into the 2030s instead of stopping at the last historical global-liquidity point.
$historical_wave_dates = $adj_dates;
$dates = dlq_monthly_dates($adj_dates[0] ?? '2000-01-01', $debt_liq_future_end);

// Preserve historical Liquidity Wave values, then carry the latest known primary-liquidity state forward.
// The forward Debt-Liquidity wave supplies the projected impulse beyond the live historical data.
$liq_wave_map = [];
for($i=0; $i<count($historical_wave_dates); $i++){
    if(isset($liq_wave_final[$i])) $liq_wave_map[$historical_wave_dates[$i]] = $liq_wave_final[$i];
}
$liq_wave_final = carry_forward($dates, $liq_wave_map);

// Re-align MOVE levels to the extended monthly timeline while retaining historical mean/stddev.
$move_vals = carry_forward($dates, $move_raw);

$ycc_raw = fetch_txt("https://crons.catenacap.xyz/experiments/Yield_Inversion_inLiquidityERA/ycc.txt");
$qe_raw  = fetch_txt("https://crons.catenacap.xyz/experiments/Yield_Inversion_inLiquidityERA/qe-tbill-issuance.txt");
$china_credit_raw = load_monthly_db_series_map($mysqli, $CHINA_CREDIT_IMPULSE_TABLE);
$ycc = carry_forward($dates, $ycc_raw);
$qe  = carry_forward($dates, $qe_raw);
$china_credit = carry_forward($dates, $china_credit_raw);
$qe_hist_last_date = dlq_last_date_in_map($qe_raw);

// ==================== LOAD GLOBAL RISK-LOVE COMPOSITE ====================
$risklove_raw = [];
$r = $mysqli->query("SELECT date, value FROM `$RISKLOVE_TABLE` ORDER BY date ASC");
if($r){
    while($row = $r->fetch_assoc()) {
        if($row['date'] === null || $row['value'] === null || !is_numeric($row['value'])) continue;
        $risklove_raw[$row['date']] = (float)$row['value'];
    }
}
$risklove_hist_last_date = dlq_last_date_in_map($risklove_raw);
$risklove = dlq_map_no_carry($dates, $risklove_raw);


// ==================== BUILD FORWARD DEBT-LIQUIDITY SERIES ====================
$debt_liq_monthly_dates = dlq_monthly_dates($dates[0] ?? '2000-01-01', $debt_liq_future_end);

$walcl_tn_raw   = dlq_fred_series_millions_to_tn('WALCL',   $FRED_API_KEY, '2000-01-01', date('Y-m-d'), $DEBT_LIQ_CACHE_DIR, $DEBT_LIQ_CACHE_TTL);
$gfdebtn_tn_raw = dlq_fred_series_millions_to_tn('GFDEBTN', $FRED_API_KEY, '2000-01-01', date('Y-m-d'), $DEBT_LIQ_CACHE_DIR, $DEBT_LIQ_CACHE_TTL);
$walcl_m        = dlq_carry_map($debt_liq_monthly_dates, $walcl_tn_raw);
$gfdebtn_m      = dlq_carry_map($debt_liq_monthly_dates, $gfdebtn_tn_raw);

$dlq_issuance_start = date('Y-m-d', strtotime('-2 years'));
$dlq_bills = dlq_fetch_db_sum($mysqli, 'US-Treasury-Security-Issuance-Bills', $dlq_issuance_start);
$dlq_notes = dlq_fetch_db_sum($mysqli, 'US-Treasury-Security-Issuance-Notes', $dlq_issuance_start);
$dlq_bonds = dlq_fetch_db_sum($mysqli, 'US-Treasury-Security-Issuance-Bonds', $dlq_issuance_start);
$dlq_rollover_weights = dlq_dynamic_rollover_weights($dlq_bills, $dlq_notes, $dlq_bonds);

$dlq_current_year = (int)date('Y');
$dlq_end_year     = $dlq_current_year + $DEBT_LIQ_FORWARD_YEARS;
$dlq_actual_maturities = dlq_fetch_treasury_maturities_by_year(2000, $dlq_end_year, $DEBT_LIQ_CACHE_DIR, $DEBT_LIQ_CACHE_TTL);
$dlq_projected_rollovers_annual_tn = dlq_project_rollovers_annual_tn($dlq_actual_maturities, $dlq_current_year, $dlq_end_year, $dlq_rollover_weights, $DEBT_LIQ_ROLLOVER_CYCLES);

$dlq_interest_tn_raw = [];
$r_dlq = $mysqli->query("SELECT record_date, interest FROM `R-Fig2-Forecasted-increase-in-Borrowing-Costs` WHERE record_date >= '2012-01-01' ORDER BY record_date ASC");
if($r_dlq){
    while($row=$r_dlq->fetch_assoc()){
        if(isset($row['record_date'], $row['interest']) && is_numeric($row['interest'])) $dlq_interest_tn_raw[$row['record_date']] = ((float)$row['interest']) / 1e12;
    }
}
$dlq_interest_m = dlq_carry_map($debt_liq_monthly_dates, $dlq_interest_tn_raw);

$debt_liq_req_stock = [];
$debt_liq_gap_tn = [];
$debt_liq_treasury_qe_need_tn = [];
$debt_liq_fed_impulse_tn = [];
$debt_liq_gross_need_tn = [];
$debt_liq_projected_debt_tn = [];
$debt_liq_stock_for_wave = [];
$dlq_last_debt = dlq_latest_value($gfdebtn_m);
$dlq_cum_interest_debt_add = 0.0;

for($i=0; $i<count($debt_liq_monthly_dates); $i++){
    $d = $debt_liq_monthly_dates[$i];
    $y = (int)substr($d,0,4);
    $histDebt = $gfdebtn_m[$d] ?? null;
    if($histDebt !== null && $d <= date('Y-m-d')){
        $projDebt = $histDebt;
        $dlq_last_debt = $histDebt;
    } else {
        $monthlyInterestAdd = (($dlq_interest_m[$d] ?? 0) / 12.0);
        $dlq_cum_interest_debt_add += max(0, $monthlyInterestAdd);
        $projDebt = ($dlq_last_debt !== null ? $dlq_last_debt : 0) + $dlq_cum_interest_debt_add;
    }

    $annualRollover = $dlq_projected_rollovers_annual_tn[$y] ?? 0.0;
    $monthlyRollover = $annualRollover / 12.0;
    $monthlyInterest = max(0, ($dlq_interest_m[$d] ?? 0) / 12.0);
    $grossNeed = $monthlyRollover + $monthlyInterest;

    $fedImpulse = null;
    if($i >= 12){
        $currWalcl = $walcl_m[$d] ?? null;
        $prevWalcl = $walcl_m[$debt_liq_monthly_dates[$i-12]] ?? null;
        $fedImpulse = ($currWalcl !== null && $prevWalcl !== null) ? max(0, ($currWalcl - $prevWalcl) / 12.0) : 0.0;
    } else {
        $fedImpulse = 0.0;
    }

    $treasuryNeed = max(0, $grossNeed - $fedImpulse);
    $reqStock = ($projDebt !== null) ? ($projDebt * $DEBT_LIQ_REQUIRED_MULTIPLE) : null;
    $walclStock = $walcl_m[$d] ?? 0.0;
    $stockForWave = $walclStock + $treasuryNeed;
    $gap = ($reqStock !== null) ? ($reqStock - $stockForWave) : null;

    $debt_liq_projected_debt_tn[] = $projDebt;
    $debt_liq_req_stock[] = $reqStock;
    $debt_liq_gap_tn[] = $gap;
    $debt_liq_treasury_qe_need_tn[] = $treasuryNeed;
    $debt_liq_fed_impulse_tn[] = $fedImpulse;
    $debt_liq_gross_need_tn[] = $grossNeed;
    $debt_liq_stock_for_wave[] = ($stockForWave > 0 ? $stockForWave : null);
}

$debt_liq_impulse = [];
for($i=12; $i<count($debt_liq_stock_for_wave); $i++){
    $a = $debt_liq_stock_for_wave[$i] ?? null;
    $b = $debt_liq_stock_for_wave[$i-12] ?? null;
    $debt_liq_impulse[] = ($a !== null && $b !== null && $a > 0 && $b > 0) ? log($a / $b) : null;
}
$debt_liq_wave_full = normalize(ema(ema(ema($debt_liq_impulse, 4), 6), 4));
$debt_liq_wave_dates = array_slice($debt_liq_monthly_dates, 12);
$debt_liq_wave_for_dates = carry_forward($dates, array_combine($debt_liq_wave_dates, $debt_liq_wave_full));

$debt_liq_plot_dates = [];
$debt_liq_wave_plot = [];
$debt_liq_gap_plot = [];
$debt_liq_treasury_need_plot = [];
$debt_liq_gross_need_plot = [];
$debt_liq_required_stock_plot = [];
$debt_liq_projected_debt_plot = [];
foreach($debt_liq_monthly_dates as $i=>$d){
    if($d < '2014-09-01') continue;
    $debt_liq_plot_dates[] = $d;
    $debt_liq_gap_plot[] = $debt_liq_gap_tn[$i] ?? null;
    $debt_liq_treasury_need_plot[] = $debt_liq_treasury_qe_need_tn[$i] ?? null;
    $debt_liq_gross_need_plot[] = $debt_liq_gross_need_tn[$i] ?? null;
    $debt_liq_required_stock_plot[] = $debt_liq_req_stock[$i] ?? null;
    $debt_liq_projected_debt_plot[] = $debt_liq_projected_debt_tn[$i] ?? null;
}
$debt_liq_wave_plot_map = array_combine($debt_liq_wave_dates, signed_log_series($debt_liq_wave_full));
foreach($debt_liq_plot_dates as $d){ $debt_liq_wave_plot[] = $debt_liq_wave_plot_map[$d] ?? null; }

$debt_liq_latest_required_stock = dlq_latest_value($debt_liq_required_stock_plot);
$debt_liq_latest_gap = dlq_latest_value($debt_liq_gap_plot);
$debt_liq_latest_treasury_need = dlq_latest_value($debt_liq_treasury_need_plot);
$debt_liq_latest_gross_need = dlq_latest_value($debt_liq_gross_need_plot);

// ==================== SECONDARY WAVES ====================
$ycc_imp=$qe_imp=$china_credit_imp=$risklove_imp=$debt_liq_imp=[];
for($i=1;$i<count($dates);$i++){
    $ycc_imp[]      = ($ycc[$i]!==null      && $ycc[$i-1]!==null)      ? ($ycc[$i]-$ycc[$i-1]) : null;
    $qe_imp[]       = ($qe[$i] !==null      && $qe[$i-1] !==null)      ? ($qe[$i] -$qe[$i-1]) : null;
    $china_credit_imp[] = ($china_credit[$i] !== null && $china_credit[$i-1] !== null) ? ($china_credit[$i] - $china_credit[$i-1]) : null;
    $risklove_imp[] = ($risklove[$i]!==null && $risklove[$i-1]!==null) ? ($risklove[$i]-$risklove[$i-1]) : null;
    $debt_liq_imp[] = ($debt_liq_wave_for_dates[$i]!==null && $debt_liq_wave_for_dates[$i-1]!==null) ? ($debt_liq_wave_for_dates[$i]-$debt_liq_wave_for_dates[$i-1]) : null;
}

$ycc_z      = normalize(ema(ema(ema($ycc_imp,4),6),4));
$qe_z       = normalize(ema(ema(ema($qe_imp ,4),6),4));
$china_credit_z = normalize(ema(ema(ema($china_credit_imp,4),6),4));
$risklove_z = normalize(ema(ema(ema($risklove_imp,4),6),4));
$debt_liq_z = normalize(ema(ema(ema($debt_liq_imp,4),6),4));

// Forward Stealth QE is projected from the debt-liquidity requirement only after the
// last actual Stealth QE source point. Global Risk-Love is deliberately not projected.
if($qe_hist_last_date !== null){
    for($j=0; $j<count($qe_z); $j++){
        $waveDate = $dates[$j+1] ?? null;
        if($waveDate !== null && $waveDate > $qe_hist_last_date && isset($debt_liq_z[$j]) && $debt_liq_z[$j] !== null){
            $qe_z[$j] = $debt_liq_z[$j];
        }
    }
}

$move_z_full = [];
foreach($move_vals as $v){
    if($v===null || $ms<=0){
        $move_z_full[] = null;
        continue;
    }
    $move_z_full[] = ($v - $mm) / $ms;
}
$move_z_full = ema($move_z_full, 3);

$risklove_move = [];
$rlCount = min(count($risklove_z), count($move_z_full));
for($i=0; $i<$rlCount; $i++){
    if($risklove_z[$i] === null){
        $risklove_move[] = null;
        continue;
    }
    $movePenalty = ($move_z_full[$i] !== null) ? (0.75 * $move_z_full[$i]) : 0.0;
    $risklove_move[] = -($risklove_z[$i] - $movePenalty);
}
$risklove_move = ema($risklove_move, 3);
$risklove_move = normalize($risklove_move);

// ==================== BTC-RELEVANT LAG ALIGNMENT ====================
$liq_btc      = shift_array($liq_wave_final, +33);
$qe_btc       = shift_array($qe_z,            0);
$ycc_btc      = shift_array($ycc_z,           0);
$china_credit_btc = shift_array($china_credit_z, $CHINA_CREDIT_IMPULSE_SHIFT_MONTHS); // China Credit Impulse impact shift is controlled in the heading settings section.
$risklove_btc = shift_array($risklove_move,   0);
$debt_liq_btc = shift_array($debt_liq_z, 0);

// ==================== DYNAMIC DOMINANCE WEIGHTS ====================
$WIN = 9;
function dominance($wave, $win){
    $d = [null];
    for($i=1;$i<count($wave);$i++){
        $d[] = ($wave[$i]!==null && $wave[$i-1]!==null)
            ? abs($wave[$i] - $wave[$i-1])
            : null;
    }
    return ema($d, $win);
}
$dom_liq      = dominance($liq_btc,      $WIN);
$dom_qe       = dominance($qe_btc ,     $WIN);
$dom_ycc      = dominance($ycc_btc,     $WIN);
$dom_china_credit = dominance($china_credit_btc, $WIN);
$dom_risklove = dominance($risklove_btc,$WIN);
$dom_debt_liq = dominance($debt_liq_btc, $WIN);

// ==================== DYNAMIC BLENDED LIQUIDITY WAVE ====================
$blend_wave = [];
$max = min(
    count($liq_wave_final),
    count($dom_liq),
    count($dom_qe),
    count($dom_ycc),
    count($dom_risklove),
    count($dom_debt_liq),
    count($dom_china_credit),
    count($liq_btc),
    count($qe_btc),
    count($ycc_btc),
    count($china_credit_btc),
    count($risklove_btc),
    count($debt_liq_btc)
);

for($i=0; $i<$max; $i++){
    if (
        !isset($dom_liq[$i], $dom_qe[$i], $dom_ycc[$i], $dom_debt_liq[$i], $dom_china_credit[$i]) ||
        !isset($liq_btc[$i], $qe_btc[$i], $ycc_btc[$i], $debt_liq_btc[$i], $china_credit_btc[$i])
    ){
        $blend_wave[] = null;
        continue;
    }

    $a = $dom_liq[$i];
    $b = $dom_qe[$i];
    $c = $dom_ycc[$i];
    $d = $dom_risklove[$i] ?? null;
    $e = $dom_debt_liq[$i];
    $f = $dom_china_credit[$i];

    if($a===null || $b===null || $c===null || $e===null || $f===null){
        $blend_wave[] = null;
        continue;
    }

    $riskValue = ($risklove_btc[$i] ?? null);
    $riskAdjDominance = ($d !== null && $riskValue !== null) ? (0.50 * $d) : 0.0;
    $debtLiqDominance = 0.75 * $e;
    $chinaCreditDominance = $CHINA_CREDIT_IMPULSE_DOMINANCE_MULTIPLE * $f;
    $sum = $a + $b + $c + $riskAdjDominance + $debtLiqDominance + $chinaCreditDominance;
    if($sum <= 0){
        $blend_wave[] = null;
        continue;
    }

    $wL = $a / $sum;
    $wQ = $b / $sum;
    $wY = $c / $sum;
    $wR = $riskAdjDominance / $sum;
    $wD = $debtLiqDominance / $sum;
    $wC = $chinaCreditDominance / $sum;

    $blend_wave[] =
        $wL * $liq_btc[$i] +
        $wQ * $qe_btc[$i] +
        $wY * $ycc_btc[$i] +
        $wR * ($riskValue ?? 0.0) +
        $wD * $debt_liq_btc[$i] +
        $wC * $china_credit_btc[$i];
}

// ==================== BTC ====================
$btc=[];
$r=$mysqli->query("SELECT date,value FROM `$BTC_TABLE` ORDER BY date ASC");
while($x=$r->fetch_assoc()) $btc[$x['date']] = floatval($x['value']);

$btc_m=[];
foreach($btc as $d=>$v){
    $btc_m[substr($d,0,7)] = $v;
}

$k = array_keys($btc_m);
sort($k);

$btc_dates = [];
$btc_vals  = [];
$btc_prices_for_hover = [];

for($i=12; $i<count($k); $i++){
    if($btc_m[$k[$i-12]] <= 0 || $btc_m[$k[$i]] <= 0) continue;

    $dt = new DateTime($k[$i].'-01');
    $btc_dates[] = $dt->format('Y-m-d');
    $btc_vals[] = 100 * log($btc_m[$k[$i]] / $btc_m[$k[$i-12]]);
    $btc_prices_for_hover[] = $btc_m[$k[$i]];
}

// ==================== FED GOV INTEREST PAYMENTS ====================
$interest_forecast = [];
$sql = "SELECT record_date, interest 
        FROM `R-Fig2-Forecasted-increase-in-Borrowing-Costs`
        WHERE record_date >= '2012-01-01'
        ORDER BY record_date";
$r = $mysqli->query($sql);
while($row = $r->fetch_assoc()){
    $interest_forecast[$row['record_date']] = floatval($row['interest']) / 1e12;
}

// ==================== PLOT FILTER ====================
$PLOT_START='2014-09-01';
$pd=$lw=$yw=$qw=$cw=$bw=$rw=$dw=[];
for($i=0;$i<count($dates);$i++){
    if($dates[$i]<$PLOT_START) continue;
    $pd[]=$dates[$i];
    $lw[]=$liq_wave_final[$i] ?? null;
    $yw[]=$ycc_z[$i-1]??null;
    $qw[]=$qe_z[$i-1]??null; // original behavior preserved
    $cw[]=$china_credit_z[$i-1]??null;
    $bw[]=$blend_wave[$i]??null;
    $rw[]=$risklove_btc[$i]??null;
    $dw[]=$debt_liq_btc[$i]??null;
}

$lw = signed_log_series($lw);
$yw = signed_log_series($yw);
$qw = signed_log_series($qw);
$cw = signed_log_series($cw);
$bw = signed_log_series($bw);
$rw = signed_log_series($rw);
$dw = signed_log_series($dw);

$bpd = $bpv = $btc_prices_hover_plot = [];
foreach($btc_dates as $i=>$d){
    if($d < $PLOT_START) continue;
    $bpd[] = $d;
    $bpv[] = $btc_vals[$i];
    $btc_prices_hover_plot[] = $btc_prices_for_hover[$i];
}


// ==================== ADD: EQUITY INDEX 12M CHANGE SERIES ====================
// Added as legend-disabled traces only. Existing BTC / macro wave calculations remain untouched.
function load_monthly_12m_log_change_from_table($mysqli, $table, $plotStart){
    $monthly = [];
    $tableEsc = $mysqli->real_escape_string($table);
    $sql = "SELECT date, value FROM `" . $tableEsc . "` WHERE date IS NOT NULL AND value IS NOT NULL ORDER BY date ASC";
    $res = $mysqli->query($sql);
    if($res){
        while($row = $res->fetch_assoc()){
            if(!isset($row['date'], $row['value']) || !is_numeric($row['value'])) continue;
            $d = substr((string)$row['date'], 0, 10);
            $monthKey = substr($d, 0, 7);
            $monthly[$monthKey] = (float)$row['value'];
        }
    }

    $keys = array_keys($monthly);
    sort($keys);

    $dates = [];
    $chg = [];
    $levels = [];
    for($i=12; $i<count($keys); $i++){
        $prev = $monthly[$keys[$i-12]] ?? null;
        $curr = $monthly[$keys[$i]] ?? null;
        if($prev === null || $curr === null || $prev <= 0 || $curr <= 0) continue;
        $dt = new DateTime($keys[$i] . '-01');
        $plotDate = $dt->format('Y-m-d');
        if($plotDate < $plotStart) continue;
        $dates[] = $plotDate;
        $chg[] = 100 * log($curr / $prev);
        $levels[] = $curr;
    }

    return [$dates, $chg, $levels];
}

[$sp500_12m_dates, $sp500_12m_vals, $sp500_12m_levels] = load_monthly_12m_log_change_from_table($mysqli, 'SP500', $PLOT_START);
[$nasdaq_12m_dates, $nasdaq_12m_vals, $nasdaq_12m_levels] = load_monthly_12m_log_change_from_table($mysqli, 'nasdaq_usd', $PLOT_START);
[$russell_12m_dates, $russell_12m_vals, $russell_12m_levels] = load_monthly_12m_log_change_from_table($mysqli, 'russell_2000_usd', $PLOT_START);
$russell_12m_dates = array_map(fn($d) => add_months_ymd($d, 5), $russell_12m_dates);

// ==================== SIGMA ENVELOPES ====================
function sigma_envelopes($wave, $window = 18){
    $d=[null];
    for($i=1;$i<count($wave);$i++) $d[] = ($wave[$i]!==null && $wave[$i-1]!==null) ? ($wave[$i]-$wave[$i-1]) : null;
    $o=['u1'=>[],'l1'=>[],'u2'=>[],'l2'=>[],'u3'=>[],'l3'=>[],'u4'=>[],'l4'=>[]];
    for($i=0;$i<count($wave);$i++){
        if($i<$window || $d[$i]===null){
            foreach($o as $kk=>$_) $o[$kk][]=null;
            continue;
        }
        $slice = array_slice($d,$i-$window+1,$window);
        $slice = array_values(array_filter($slice, fn($v)=>$v!==null));
        $s=stddev($slice);
        foreach([1,2,3,4] as $n){
            $o["u$n"][] = ($wave[$i]!==null) ? ($wave[$i]+$n*$s) : null;
            $o["l$n"][] = ($wave[$i]!==null) ? ($wave[$i]-$n*$s) : null;
        }
    }
    return $o;
}
$liq_sig      = sigma_envelopes($lw);
$ycc_sig      = sigma_envelopes($yw);
$qe_sig       = sigma_envelopes($qw);
$china_credit_sig = sigma_envelopes($cw);
$china_credit_plot_dates = array_map(fn($d) => add_months_ymd($d, $CHINA_CREDIT_IMPULSE_SHIFT_MONTHS), $pd);
$blend_sig    = sigma_envelopes($bw);
$risklove_sig = sigma_envelopes($rw);
$debt_liq_sig = sigma_envelopes($dw);

// ==================== ADD: MOVE RAW DAILY ====================
$move_raw_dates=[];
$move_raw_vals=[];
foreach($move_raw as $d=>$v){
    if ($d < $PLOT_START) continue;
    $move_raw_dates[] = $d;
    $move_raw_vals[]  = $v;
}

// ==================== VISUAL ANALYSIS SETTINGS ====================
$VISUAL_COMPARE_FROM = '2019-08-01';
$SHIFT_BTC_PLOT      =  2;
$SHIFT_QE_PLOT       = 10;   // same original visual shift
$SHIFT_BLEND_PLOT    = 10;
$SHIFT_RISKLOVE_PLOT = 10;

$btc_plot_map      = build_shifted_map($btc_dates, $btc_vals, $SHIFT_BTC_PLOT,      $VISUAL_COMPARE_FROM);
$qe_plot_map       = build_shifted_map($pd ?? [], $qw ?? [], $SHIFT_QE_PLOT,        $VISUAL_COMPARE_FROM);
$blend_plot_map    = build_shifted_map($pd ?? [], $bw ?? [], $SHIFT_BLEND_PLOT,     $VISUAL_COMPARE_FROM);
$risklove_plot_map = build_shifted_map($pd ?? [], $rw ?? [], $SHIFT_RISKLOVE_PLOT,  $VISUAL_COMPARE_FROM);

// ==================== IMPACT ANALYSIS ====================
[$dates_qe,        $btc_vs_qe_btc,        $btc_vs_qe_wave]        = align_two_series($btc_plot_map, $qe_plot_map,       $VISUAL_COMPARE_FROM);
[$dates_blend,     $btc_vs_blend_btc,     $btc_vs_blend_wave]     = align_two_series($btc_plot_map, $blend_plot_map,    $VISUAL_COMPARE_FROM);
[$dates_risklove,  $btc_vs_risklove_btc,  $btc_vs_risklove_wave]  = align_two_series($btc_plot_map, $risklove_plot_map, $VISUAL_COMPARE_FROM);

$qe_corr       = pearson_corr($btc_vs_qe_btc,       $btc_vs_qe_wave);
$blend_corr    = pearson_corr($btc_vs_blend_btc,    $btc_vs_blend_wave);
$risklove_corr = pearson_corr($btc_vs_risklove_btc, $btc_vs_risklove_wave);

$qe_r2         = r_squared_from_corr($qe_corr);
$blend_r2      = r_squared_from_corr($blend_corr);
$risklove_r2   = r_squared_from_corr($risklove_corr);

$qe_dir_hit       = directional_hit_rate($btc_vs_qe_btc,       $btc_vs_qe_wave);
$blend_dir_hit    = directional_hit_rate($btc_vs_blend_btc,    $btc_vs_blend_wave);
$risklove_dir_hit = directional_hit_rate($btc_vs_risklove_btc, $btc_vs_risklove_wave);

$qe_turn_hit       = turning_point_hit_rate($btc_vs_qe_btc,       $btc_vs_qe_wave);
$blend_turn_hit    = turning_point_hit_rate($btc_vs_blend_btc,    $btc_vs_blend_wave);
$risklove_turn_hit = turning_point_hit_rate($btc_vs_risklove_btc, $btc_vs_risklove_wave);

$qe_impact_score       = composite_impact_score($qe_corr,       $qe_r2,       $qe_dir_hit,       $qe_turn_hit);
$blend_impact_score    = composite_impact_score($blend_corr,    $blend_r2,    $blend_dir_hit,    $blend_turn_hit);
$risklove_impact_score = composite_impact_score($risklove_corr, $risklove_r2, $risklove_dir_hit, $risklove_turn_hit);

$impactScores = [
    'Blended Wave' => $blend_impact_score,
    'Stealth QE Wave' => $qe_impact_score,
    'Global Risk-Love Wave' => $risklove_impact_score
];
arsort($impactScores);
$impact_winner = array_key_first($impactScores);

// ==================== ROLLING DOMINANCE REGIME ====================
$common_all_dates = array_values(array_intersect(array_keys($btc_plot_map), array_keys($qe_plot_map), array_keys($blend_plot_map), array_keys($risklove_plot_map)));
sort($common_all_dates);

$btc_all      = [];
$qe_all       = [];
$blend_all    = [];
$risklove_all = [];
foreach($common_all_dates as $d){
    $btc_all[]      = $btc_plot_map[$d];
    $qe_all[]       = $qe_plot_map[$d];
    $blend_all[]    = $blend_plot_map[$d];
    $risklove_all[] = $risklove_plot_map[$d];
}

$rollWin = 12;
$rollingCorrQE = [];
$rollingCorrBlend = [];
$rollingCorrRiskLove = [];
$dominanceState = [];

for($i=0; $i<count($common_all_dates); $i++){
    if($i < $rollWin - 1){
        $rollingCorrQE[] = null;
        $rollingCorrBlend[] = null;
        $rollingCorrRiskLove[] = null;
        $dominanceState[] = 'neutral';
        continue;
    }

    $btcSlice      = array_slice($btc_all,      $i-$rollWin+1, $rollWin);
    $qeSlice       = array_slice($qe_all,       $i-$rollWin+1, $rollWin);
    $blendSlice    = array_slice($blend_all,    $i-$rollWin+1, $rollWin);
    $riskloveSlice = array_slice($risklove_all, $i-$rollWin+1, $rollWin);

    $cq = pearson_corr($btcSlice, $qeSlice);
    $cb = pearson_corr($btcSlice, $blendSlice);
    $cr = pearson_corr($btcSlice, $riskloveSlice);

    $rollingCorrQE[] = $cq;
    $rollingCorrBlend[] = $cb;
    $rollingCorrRiskLove[] = $cr;

    $scores = [];
    if($cb !== null) $scores['blend'] = $cb;
    if($cq !== null) $scores['qe'] = $cq;
    if($cr !== null) $scores['risklove'] = $cr;

    if(count($scores) >= 2){
        arsort($scores);
        $keys = array_keys($scores);
        $vals = array_values($scores);
        $dominanceState[] = (($vals[0] - $vals[1]) > 0.05) ? $keys[0] : 'neutral';
    } else {
        $dominanceState[] = 'neutral';
    }
}
$dominance_regions = compress_regimes($common_all_dates, $dominanceState);

// ==================== COMMON ARRAYS ====================
$commonDates3 = $common_all_dates;
$btc3      = $btc_all;
$qe3       = $qe_all;
$blend3    = $blend_all;
$risklove3 = $risklove_all;

// ==================== PEAK ANALYSIS ====================
$btc_cycle_peaks_actual = find_btc_cycle_peaks($btc_m, 5, 5, 12, 0.25, 8);
$btc_cycle_peaks_plot = [];
foreach($btc_cycle_peaks_actual as $d){
    $btc_cycle_peaks_plot[] = add_months_ymd($d . '-01', $SHIFT_BTC_PLOT);
}

$peak_rows = [];
$peakWinsQE = 0;
$peakWinsBlend = 0;
$qePeakScores = [];
$blendPeakScores = [];

$qePeakConditionCounts = [
    'btc_hot'=>0,'wave_rollover_3m'=>0,'wave_below_3m_high'=>0,'wave_below_6m_ma'=>0,'btc_up_wave_down_2m'=>0,'corr6_lt_corr12'=>0
];
$blendPeakConditionCounts = [
    'btc_hot'=>0,'wave_rollover_3m'=>0,'wave_below_3m_high'=>0,'wave_below_6m_ma'=>0,'btc_up_wave_down_2m'=>0,'corr6_lt_corr12'=>0
];
$peakObsCount = 0;

foreach($btc_cycle_peaks_plot as $peakPlotDate){
    $idx = array_search($peakPlotDate, $commonDates3, true);
    if($idx === false || $idx < 12) continue;

    $qeRes    = score_peak_driver_for_wave($peakPlotDate, $commonDates3, $btc3, $qe3);
    $blendRes = score_peak_driver_for_wave($peakPlotDate, $commonDates3, $btc3, $blend3);

    if(!$qeRes || !$blendRes) continue;

    $winner = ($blendRes['peakScore'] > $qeRes['peakScore']) ? 'Blended Wave' : 'Stealth QE Wave';

    if($winner === 'Blended Wave') $peakWinsBlend++;
    else $peakWinsQE++;

    $qePeakScores[]    = $qeRes['peakScore'];
    $blendPeakScores[] = $blendRes['peakScore'];

    $winStart = max(0, $idx - 11);
    $btcWin   = array_slice($btc3,   $winStart, 12);
    $qeWin    = array_slice($qe3,    $winStart, 12);
    $blendWin = array_slice($blend3, $winStart, 12);
    $dWin     = array_slice($commonDates3, $winStart, 12);

    $qeCond    = condition_pack_peak($dWin, $btcWin, $qeWin);
    $blendCond = condition_pack_peak($dWin, $btcWin, $blendWin);

    foreach($qePeakConditionCounts as $kCond => $v){
        if($qeCond[$kCond]) $qePeakConditionCounts[$kCond]++;
    }
    foreach($blendPeakConditionCounts as $kCond => $v){
        if($blendCond[$kCond]) $blendPeakConditionCounts[$kCond]++;
    }

    $peakObsCount++;

    $peak_rows[] = [
        'date'            => $peakPlotDate,
        'qeScore'         => $qeRes['peakScore'],
        'blendScore'      => $blendRes['peakScore'],
        'qeLeadMonths'    => $qeRes['leadMonths'],
        'blendLeadMonths' => $blendRes['leadMonths'],
        'winner'          => $winner
    ];
}

$avgQePeakScore    = count($qePeakScores)    ? mean($qePeakScores)    : null;
$avgBlendPeakScore = count($blendPeakScores) ? mean($blendPeakScores) : null;
$peak_driver_winner = ($avgBlendPeakScore > $avgQePeakScore) ? 'Blended Wave' : 'Stealth QE Wave';

$currentQePeakCond    = condition_pack_peak($commonDates3, $btc3, $qe3);
$currentBlendPeakCond = condition_pack_peak($commonDates3, $btc3, $blend3);

$qePeakTemplateRates = [];
$blendPeakTemplateRates = [];
foreach($qePeakConditionCounts as $kk => $v){
    $qePeakTemplateRates[$kk] = $peakObsCount ? ($v / $peakObsCount) : 0;
}
foreach($blendPeakConditionCounts as $kk => $v){
    $blendPeakTemplateRates[$kk] = $peakObsCount ? ($v / $peakObsCount) : 0;
}
$qePeakTemplateEval    = evaluate_template_match($currentQePeakCond, $qePeakTemplateRates, 0.60);
$blendPeakTemplateEval = evaluate_template_match($currentBlendPeakCond, $blendPeakTemplateRates, 0.60);

$current_peak_condition_winner = ($blendPeakTemplateEval['met'] > $qePeakTemplateEval['met']) ? 'Blended Wave'
    : (($qePeakTemplateEval['met'] > $blendPeakTemplateEval['met']) ? 'Stealth QE Wave' : 'Tie');

// ==================== BOTTOM ANALYSIS ====================
$btc_cycle_bottoms_actual = find_btc_cycle_bottoms($btc_m, 5, 5, 12, 0.35, 8);
$btc_cycle_bottoms_plot = [];
foreach($btc_cycle_bottoms_actual as $d){
    $btc_cycle_bottoms_plot[] = add_months_ymd($d . '-01', $SHIFT_BTC_PLOT);
}

$bottom_rows = [];
$bottomWinsQE = 0;
$bottomWinsBlend = 0;
$qeBottomScores = [];
$blendBottomScores = [];

$qeBottomConditionCounts = [
    'btc_cold'=>0,'wave_turnup_3m'=>0,'wave_above_3m_low'=>0,'wave_above_6m_ma'=>0,'btc_down_wave_up_2m'=>0,'corr6_gt_corr12'=>0
];
$blendBottomConditionCounts = [
    'btc_cold'=>0,'wave_turnup_3m'=>0,'wave_above_3m_low'=>0,'wave_above_6m_ma'=>0,'btc_down_wave_up_2m'=>0,'corr6_gt_corr12'=>0
];
$bottomObsCount = 0;

foreach($btc_cycle_bottoms_plot as $bottomPlotDate){
    $idx = array_search($bottomPlotDate, $commonDates3, true);
    if($idx === false || $idx < 12) continue;

    $qeRes    = score_bottom_driver_for_wave($bottomPlotDate, $commonDates3, $btc3, $qe3);
    $blendRes = score_bottom_driver_for_wave($bottomPlotDate, $commonDates3, $btc3, $blend3);

    if(!$qeRes || !$blendRes) continue;

    $winner = ($blendRes['bottomScore'] > $qeRes['bottomScore']) ? 'Blended Wave' : 'Stealth QE Wave';

    if($winner === 'Blended Wave') $bottomWinsBlend++;
    else $bottomWinsQE++;

    $qeBottomScores[]    = $qeRes['bottomScore'];
    $blendBottomScores[] = $blendRes['bottomScore'];

    $winStart = max(0, $idx - 11);
    $btcWin   = array_slice($btc3,   $winStart, 12);
    $qeWin    = array_slice($qe3,    $winStart, 12);
    $blendWin = array_slice($blend3, $winStart, 12);
    $dWin     = array_slice($commonDates3, $winStart, 12);

    $qeCond    = condition_pack_bottom($dWin, $btcWin, $qeWin);
    $blendCond = condition_pack_bottom($dWin, $btcWin, $blendWin);

    foreach($qeBottomConditionCounts as $kCond => $v){
        if($qeCond[$kCond]) $qeBottomConditionCounts[$kCond]++;
    }
    foreach($blendBottomConditionCounts as $kCond => $v){
        if($blendCond[$kCond]) $blendBottomConditionCounts[$kCond]++;
    }

    $bottomObsCount++;

    $bottom_rows[] = [
        'date'            => $bottomPlotDate,
        'qeScore'         => $qeRes['bottomScore'],
        'blendScore'      => $blendRes['bottomScore'],
        'qeLeadMonths'    => $qeRes['leadMonths'],
        'blendLeadMonths' => $blendRes['leadMonths'],
        'winner'          => $winner
    ];
}

$avgQeBottomScore    = count($qeBottomScores)    ? mean($qeBottomScores)    : null;
$avgBlendBottomScore = count($blendBottomScores) ? mean($blendBottomScores) : null;
$bottom_driver_winner = ($avgBlendBottomScore > $avgQeBottomScore) ? 'Blended Wave' : 'Stealth QE Wave';

$currentQeBottomCond    = condition_pack_bottom($commonDates3, $btc3, $qe3);
$currentBlendBottomCond = condition_pack_bottom($commonDates3, $btc3, $blend3);

$qeBottomTemplateRates = [];
$blendBottomTemplateRates = [];
foreach($qeBottomConditionCounts as $kk => $v){
    $qeBottomTemplateRates[$kk] = $bottomObsCount ? ($v / $bottomObsCount) : 0;
}
foreach($blendBottomConditionCounts as $kk => $v){
    $blendBottomTemplateRates[$kk] = $bottomObsCount ? ($v / $bottomObsCount) : 0;
}
$qeBottomTemplateEval    = evaluate_template_match($currentQeBottomCond, $qeBottomTemplateRates, 0.60);
$blendBottomTemplateEval = evaluate_template_match($currentBlendBottomCond, $blendBottomTemplateRates, 0.60);

$current_bottom_condition_winner = ($blendBottomTemplateEval['met'] > $qeBottomTemplateEval['met']) ? 'Blended Wave'
    : (($qeBottomTemplateEval['met'] > $blendBottomTemplateEval['met']) ? 'Stealth QE Wave' : 'Tie');

// ==================== CURRENT ROLLING CORRS ====================
$latestRollingQE       = null;
$latestRollingBlend    = null;
$latestRollingRiskLove = null;
for($i=count($rollingCorrQE)-1; $i>=0; $i--){
    if($rollingCorrQE[$i] !== null){ $latestRollingQE = $rollingCorrQE[$i]; break; }
}
for($i=count($rollingCorrBlend)-1; $i>=0; $i--){
    if($rollingCorrBlend[$i] !== null){ $latestRollingBlend = $rollingCorrBlend[$i]; break; }
}
for($i=count($rollingCorrRiskLove)-1; $i>=0; $i--){
    if($rollingCorrRiskLove[$i] !== null){ $latestRollingRiskLove = $rollingCorrRiskLove[$i]; break; }
}

// ==================== MARKERS ====================
$btc_peak_marker_dates = [];
$btc_peak_marker_vals  = [];
foreach($btc_cycle_peaks_plot as $pdPeak){
    if(isset($btc_plot_map[$pdPeak])){
        $btc_peak_marker_dates[] = $pdPeak;
        $btc_peak_marker_vals[]  = $btc_plot_map[$pdPeak];
    }
}

$btc_bottom_marker_dates = [];
$btc_bottom_marker_vals  = [];
foreach($btc_cycle_bottoms_plot as $pdBottom){
    if(isset($btc_plot_map[$pdBottom])){
        $btc_bottom_marker_dates[] = $pdBottom;
        $btc_bottom_marker_vals[]  = $btc_plot_map[$pdBottom];
    }
}

?>
<!DOCTYPE html>
<html lang="en" data-bs-theme="dark">
<head>
  <meta charset="UTF-8">
  <meta name="viewport" content="width=device-width, initial-scale=1.0">

  <link href="https://fonts.googleapis.com/css?family=Nunito" rel="stylesheet" />
  <link rel="stylesheet" href="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/css/bootstrap.min.css" crossorigin="anonymous">
  <link href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.6.0/css/all.min.css" rel="stylesheet" />
  <link href="../style.css?v=<?php echo time(); ?>" rel="stylesheet" />

  <script type="text/javascript" src="https://code.jquery.com/jquery-3.7.1.min.js"></script>
  <script type="text/javascript" src="https://code.jquery.com/jquery-migrate-3.5.2.min.js"></script>
  <script type="text/javascript" src="https://code.jquery.com/ui/1.14.0/jquery-ui.min.js"></script>
  <script src="https://cdn.jsdelivr.net/npm/bootstrap@5.3.2/dist/js/bootstrap.min.js" crossorigin="anonymous"></script>
  <script src="https://cdn.plot.ly/plotly-2.30.0.min.js"></script>

  <style>
    body{ background-color:#212529; color:#fff; }
    #chart{ width:100%; height:500px; }
    .card {font-size: 14px;}
    .table td, .table th{ font-size: 14px; white-space:nowrap; }
  </style>
</head>
<body>

	<div class="container-fluid py-3">
	  <div class="row g-3 mb-3">
	    <div class="col-lg-3 col-sm-6">
	      <div class="card text-light h-100">
	        <div class="card-body">
	          <h5 class="card-title mb-3">Impact Since Dec August 2019</h5>
	          <div><strong>Visual winner:</strong> <?=htmlspecialchars($impact_winner)?></div>
	          <div><strong>Blended impact score:</strong> <?=number_format((float)$blend_impact_score, 2)?></div>
	          <div><strong>Stealth QE impact score:</strong> <?=number_format((float)$qe_impact_score, 2)?></div>
	          <hr>
	          <div><strong>Blended corr:</strong> <?=($blend_corr!==null?number_format($blend_corr,3):'n/a')?></div>
	          <div><strong>QE corr:</strong> <?=($qe_corr!==null?number_format($qe_corr,3):'n/a')?></div>
	          <div><strong>Risk-Love corr:</strong> <?=($risklove_corr!==null?number_format($risklove_corr,3):'n/a')?></div>
	          <div><strong>Blended R²:</strong> <?=($blend_r2!==null?number_format($blend_r2,3):'n/a')?></div>
	          <div><strong>QE R²:</strong> <?=($qe_r2!==null?number_format($qe_r2,3):'n/a')?></div>
	          <div><strong>Risk-Love R²:</strong> <?=($risklove_r2!==null?number_format($risklove_r2,3):'n/a')?></div>
	        </div>
	      </div>
	    </div>
	
	    <div class="col-lg-3 col-sm-6">
	      <div class="card text-light h-100">
	        <div class="card-body">
	          <h5 class="card-title mb-3">BTC Peak Driver</h5>
	          <div><strong>Historical winner:</strong> <?=htmlspecialchars($peak_driver_winner)?></div>
	          <div><strong>Blended avg peak score:</strong> <?=($avgBlendPeakScore!==null?number_format($avgBlendPeakScore,2):'n/a')?></div>
	          <div><strong>QE avg peak score:</strong> <?=($avgQePeakScore!==null?number_format($avgQePeakScore,2):'n/a')?></div>
	          <div><strong>Blended peak wins:</strong> <?=$peakWinsBlend?></div>
	          <div><strong>QE peak wins:</strong> <?=$peakWinsQE?></div>
	          <div><strong>Peaks detected:</strong> <?=$peakObsCount?></div>
	          <hr>
	          <div><strong>Current peak template:</strong></div>
	          <div>Blended: <?=htmlspecialchars($blendPeakTemplateEval['status'])?> (<?=$blendPeakTemplateEval['met']?>/<?=$blendPeakTemplateEval['need']?>)</div>
	          <div>QE: <?=htmlspecialchars($qePeakTemplateEval['status'])?> (<?=$qePeakTemplateEval['met']?>/<?=$qePeakTemplateEval['need']?>)</div>
	          <div><strong>Current peak leader:</strong> <?=htmlspecialchars($current_peak_condition_winner)?></div>
	        </div>
	      </div>
	    </div>
	
	    <div class="col-lg-3 col-sm-6">
	      <div class="card text-light h-100">
	        <div class="card-body">
	          <h5 class="card-title mb-3">BTC Bottom Driver</h5>
	          <div><strong>Historical winner:</strong> <?=htmlspecialchars($bottom_driver_winner)?></div>
	          <div><strong>Blended avg bottom score:</strong> <?=($avgBlendBottomScore!==null?number_format($avgBlendBottomScore,2):'n/a')?></div>
	          <div><strong>QE avg bottom score:</strong> <?=($avgQeBottomScore!==null?number_format($avgQeBottomScore,2):'n/a')?></div>
	          <div><strong>Blended bottom wins:</strong> <?=$bottomWinsBlend?></div>
	          <div><strong>QE bottom wins:</strong> <?=$bottomWinsQE?></div>
	          <div><strong>Bottoms detected:</strong> <?=$bottomObsCount?></div>
	          <hr>
	          <div><strong>Current bottom template:</strong></div>
	          <div>Blended: <?=htmlspecialchars($blendBottomTemplateEval['status'])?> (<?=$blendBottomTemplateEval['met']?>/<?=$blendBottomTemplateEval['need']?>)</div>
	          <div>QE: <?=htmlspecialchars($qeBottomTemplateEval['status'])?> (<?=$qeBottomTemplateEval['met']?>/<?=$qeBottomTemplateEval['need']?>)</div>
	          <div><strong>Current bottom leader:</strong> <?=htmlspecialchars($current_bottom_condition_winner)?></div>
	        </div>
	      </div>
	    </div>
	
	    <div class="col-lg-3 col-sm-6">
	      <div class="card text-light h-100">
	        <div class="card-body">
	          <h5 class="card-title mb-3">Current Conditions</h5>
	          <div><strong>Latest 12m rolling corr:</strong></div>
	          <div>Blended: <?=($latestRollingBlend!==null?number_format($latestRollingBlend,3):'n/a')?></div>
	          <div>QE: <?=($latestRollingQE!==null?number_format($latestRollingQE,3):'n/a')?></div>
	          <div>Risk-Love: <?=($latestRollingRiskLove!==null?number_format($latestRollingRiskLove,3):'n/a')?></div>
	          <hr>
	          <div><strong>Primary peak setup leader:</strong> <?=htmlspecialchars($current_peak_condition_winner)?></div>
	          <div><strong>Primary bottom setup leader:</strong> <?=htmlspecialchars($current_bottom_condition_winner)?></div>
	          <div><strong>Net impact winner:</strong> <?=htmlspecialchars($impact_winner)?></div>
	        </div>
	      </div>
	    </div>
	  </div>
	


      <div class="row g-3 mb-3">
        <div class="col-12">
          <div class="card text-light h-100">
            <div class="card-body">
              <h5 class="card-title mb-3">Debt-Liquidity Forward Roadmap</h5>
              <div class="row">
                <div class="col-md-3"><strong>Required liquidity multiple:</strong> <?=number_format($DEBT_LIQ_REQUIRED_MULTIPLE,2)?>x debt</div>
                <div class="col-md-3"><strong>Latest required stock:</strong> <?=($debt_liq_latest_required_stock!==null?'$'.number_format($debt_liq_latest_required_stock,2).'T':'n/a')?></div>
                <div class="col-md-3"><strong>Latest liquidity gap:</strong> <?=($debt_liq_latest_gap!==null?'$'.number_format($debt_liq_latest_gap,2).'T':'n/a')?></div>
                <div class="col-md-3"><strong>Latest Treasury QE need:</strong> <?=($debt_liq_latest_treasury_need!==null?'$'.number_format($debt_liq_latest_treasury_need,3).'T/mo':'n/a')?></div>
              </div>
              <div class="mt-2 small text-secondary">
                Rollover projection uses current/future Treasury maturities, recent Bills/Notes/Bonds issuance mix, WALCL and GFDEBTN from FRED, and forecast interest outlay where available.
              </div>
            </div>
          </div>
        </div>
      </div>

	  <div id="chart"></div>
	
	  <?php if(!empty($peak_rows)){ ?>
	  <div class="row g-3 mt-2">
	    <div class="col-12">
	      <div class="card text-light">
	        <div class="card-body">
	          <h5 class="card-title mb-3">Historical BTC Peak Diagnostics</h5>
	          <div class="table-responsive">
	            <table class="table table-dark table-sm table-striped align-middle mb-0">
	              <thead>
	                <tr>
	                  <th>Peak (plotted)</th>
	                  <th>QE Score</th>
	                  <th>QE Lead (m)</th>
	                  <th>Blended Score</th>
	                  <th>Blended Lead (m)</th>
	                  <th>Winner</th>
	                </tr>
	              </thead>
	              <tbody>
	                <?php foreach($peak_rows as $row){ ?>
	                <tr>
	                  <td><?=htmlspecialchars($row['date'])?></td>
	                  <td><?=number_format($row['qeScore'],2)?></td>
	                  <td><?=htmlspecialchars((string)$row['qeLeadMonths'])?></td>
	                  <td><?=number_format($row['blendScore'],2)?></td>
	                  <td><?=htmlspecialchars((string)$row['blendLeadMonths'])?></td>
	                  <td><strong><?=htmlspecialchars($row['winner'])?></strong></td>
	                </tr>
	                <?php } ?>
	              </tbody>
	            </table>
	          </div>
	        </div>
	      </div>
	    </div>
	  </div>
	  <?php } ?>
	
	  <?php if(!empty($bottom_rows)){ ?>
	  <div class="row g-3 mt-3">
	    <div class="col-12">
	      <div class="card text-light">
	        <div class="card-body">
	          <h5 class="card-title mb-3">Historical BTC Bottom Diagnostics</h5>
	          <div class="table-responsive">
	            <table class="table table-dark table-sm table-striped align-middle mb-0">
	              <thead>
	                <tr>
	                  <th>Bottom (plotted)</th>
	                  <th>QE Score</th>
	                  <th>QE Lead (m)</th>
	                  <th>Blended Score</th>
	                  <th>Blended Lead (m)</th>
	                  <th>Winner</th>
	                </tr>
	              </thead>
	              <tbody>
	                <?php foreach($bottom_rows as $row){ ?>
	                <tr>
	                  <td><?=htmlspecialchars($row['date'])?></td>
	                  <td><?=number_format($row['qeScore'],2)?></td>
	                  <td><?=htmlspecialchars((string)$row['qeLeadMonths'])?></td>
	                  <td><?=number_format($row['blendScore'],2)?></td>
	                  <td><?=htmlspecialchars((string)$row['blendLeadMonths'])?></td>
	                  <td><strong><?=htmlspecialchars($row['winner'])?></strong></td>
	                </tr>
	                <?php } ?>
	              </tbody>
	            </table>
	          </div>
	        </div>
	      </div>
	    </div>
	  </div>
	  <?php } ?>
	</div>
	
	<script>
		const dominanceShapes = [
		<?php foreach($dominance_regions as $r){
		    if($r['state'] === 'neutral') continue;
		    $fill = ($r['state'] === 'blend') ? 'rgba(100,181,255,0.08)' : (($r['state'] === 'risklove') ? 'rgba(255,105,180,0.08)' : 'rgba(124,252,152,0.08)');
		?>
		{
		  type:'rect',
		  xref:'x',
		  yref:'paper',
		  x0:'<?= $r['start'] ?>',
		  x1:'<?= $r['end'] ?>',
		  y0:0,
		  y1:1,
		  line:{width:0},
		  fillcolor:'<?= $fill ?>',
		  layer:'below'
		},
		<?php } ?>
		];
		
		const traces = [
		  {
		    x: <?=json_encode($move_raw_dates)?>,
		    y: <?=json_encode($move_raw_vals)?>,
		    type:'bar',
		    yaxis:'y3',
		    name:'MOVE Index',
		    opacity:0.25
		  },
		
		  <?php foreach([4,3,2,1] as $k){ ?>
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($liq_sig["u$k"])?>,
		    legendgroup:'liq',
		    showlegend:false,
		    visible:'legendonly',
		    line:{color:'rgba(181,126,220,<?=0.025*$k?>)'}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($liq_sig["l$k"])?>,
		    legendgroup:'liq',
		    fill:'tonexty',
		    showlegend:false,
		    visible:'legendonly',
		    line:{color:'rgba(181,126,220,<?=0.025*$k?>)'}
		  },
		  <?php } ?>
		
		  <?php foreach([4,3,2,1] as $k){ ?>
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($ycc_sig["u$k"])?>,
		    legendgroup:'ycc',
		    showlegend:false,
		    visible:'legendonly',
		    line:{color:'rgba(255,179,71,<?=0.025*$k?>)'}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($ycc_sig["l$k"])?>,
		    legendgroup:'ycc',
		    fill:'tonexty',
		    showlegend:false,
		    visible:'legendonly',
		    line:{color:'rgba(255,179,71,<?=0.025*$k?>)'}
		  },
		  <?php } ?>
          /* ===== QE deviation bubble restored and visible ===== */
          <?php foreach([4,3,2,1] as $k){ ?>
          {
            x: <?=json_encode($pd)?>,
            y: <?=json_encode($qe_sig["u$k"])?>,
            legendgroup:'qe',
            showlegend:false,
            visible:true,
            mode:'lines',
            hoverinfo:'skip',
            line:{color:'rgba(124,252,152,0)', width:0.5}
          },
          {
            x: <?=json_encode($pd)?>,
            y: <?=json_encode($qe_sig["l$k"])?>,
            legendgroup:'qe',
            fill:'tonexty',
            fillcolor:'rgba(124,252,152,<?=0.035*$k?>)',
            showlegend:false,
            visible:true,
            mode:'lines',
            hoverinfo:'skip',
            line:{color:'rgba(124,252,152,0)', width:0.5}
          },
          <?php } ?>
		
          /* ===== China Credit Impulse Stealth QE-style deviation bubble, disabled by default ===== */
          <?php foreach([4,3,2,1] as $k){ ?>
          {
            x: <?=json_encode($china_credit_plot_dates)?>,
            y: <?=json_encode($china_credit_sig["u$k"])?>,
            legendgroup:'china_credit_impulse',
            showlegend:false,
            visible:'legendonly',
            mode:'lines',
            hoverinfo:'skip',
            line:{color:'rgba(245,255,0,0)', width:0.5}
          },
          {
            x: <?=json_encode($china_credit_plot_dates)?>,
            y: <?=json_encode($china_credit_sig["l$k"])?>,
            legendgroup:'china_credit_impulse',
            fill:'tonexty',
            fillcolor:'rgba(245,255,0,<?=0.04*$k?>)',
            showlegend:false,
            visible:'legendonly',
            mode:'lines',
            hoverinfo:'skip',
            line:{color:'rgba(245,255,0,0)', width:0.5}
          },
          <?php } ?>

		  <?php foreach([4,3,2,1] as $k){ ?>
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($blend_sig["u$k"])?>,
		    legendgroup:'blend',
		    showlegend:false,
		    line:{color:'rgba(100,181,255,<?=0.025*$k?>)'}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($blend_sig["l$k"])?>,
		    legendgroup:'blend',
		    fill:'tonexty',
		    showlegend:false,
		    line:{color:'rgba(100,181,255,<?=0.025*$k?>)'}
		  },
		  <?php } ?>

		  <?php foreach([4,3,2,1] as $k){ ?>
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($debt_liq_sig["u$k"])?>,
		    legendgroup:'debtliq',
		    showlegend:false,
		    visible:true,
		    line:{color:'rgba(0,229,255,<?=0.025*$k?>)'}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($debt_liq_sig["l$k"])?>,
		    legendgroup:'debtliq',
		    fill:'tonexty',
		    showlegend:false,
		    visible:true,
		    line:{color:'rgba(0,229,255,<?=0.025*$k?>)'}
		  },
		  <?php } ?>

		  <?php foreach([4,3,2,1] as $k){ ?>
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($risklove_sig["u$k"])?>,
		    legendgroup:'risklove',
		    showlegend:false,
		    visible:'legendonly',
		    line:{color:'rgba(255,105,180,<?=0.10*$k?>)'}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($risklove_sig["l$k"])?>,
		    legendgroup:'risklove',
		    fill:'tonexty',
		    showlegend:false,
		    visible:'legendonly',
		    line:{color:'rgba(255,105,180,<?=0.10*$k?>)'}
		  },
		  <?php } ?>
		
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($lw)?>,
		    name:'Liquidity Wave',
		    legendgroup:'liq',
		    visible:'legendonly',
		    line:{color:'#b57edc',width:2}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($yw)?>,
		    name:'YCC Wave',
		    legendgroup:'ycc',
		    visible:'legendonly',
		    line:{color:'#ffb347',width:1.8}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($qw)?>,
		    name:'Stealth QE Wave',
		    legendgroup:'qe',
		    visible:true,
		    line:{color:'#7CFC98',width:1.4}
		  },
          {
            x: <?=json_encode($china_credit_plot_dates)?>,
            y: <?=json_encode($cw)?>,
            name:'China Credit Impulse Stealth QE Process (+' + <?=json_encode($CHINA_CREDIT_IMPULSE_SHIFT_MONTHS)?> + 'M)',
            legendgroup:'china_credit_impulse',
            visible:'legendonly',
            line:{color:'#f5ff00',width:2.1,dash:'dash'}
          },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($bw)?>,
		    name:'Blended Wave',
		    legendgroup:'blend',
		    visible:true,
		    line:{color:'#64b5ff',width:1.9}
		  },
		  {
		    x: <?=json_encode($pd)?>,
		    y: <?=json_encode($rw)?>,
		    name:'Global Risk-Love Wave Inverted',
		    legendgroup:'risklove',
		    visible:'legendonly',
		    line:{color:'#ff69b4',width:3.4}
		  },
          {
            x: <?=json_encode($pd)?>,
            y: <?=json_encode($dw)?>,
            name:'Debt-Liquidity Requirement Wave',
            legendgroup:'debtliq',
            visible:true,
            line:{color:'#00e5ff',width:2.4,dash:'dot'}
          },
          {
            x: <?=json_encode($debt_liq_plot_dates)?>,
            y: <?=json_encode($debt_liq_treasury_need_plot)?>,
            name:'Treasury QE Need / Liquidity Gap Fill ($T/mo)',
            yaxis:'y4',
            visible:'legendonly',
            line:{color:'#00ffc8',width:1.8}
          },
          {
            x: <?=json_encode($debt_liq_plot_dates)?>,
            y: <?=json_encode($debt_liq_gross_need_plot)?>,
            name:'Gross Debt Liquidity Need ($T/mo)',
            yaxis:'y4',
            visible:'legendonly',
            line:{color:'#ffd166',width:1.4}
          },
          {
            x: <?=json_encode($debt_liq_plot_dates)?>,
            y: <?=json_encode($debt_liq_required_stock_plot)?>,
            name:'Required Liquidity Stock 270% Debt ($T)',
            yaxis:'y4',
            visible:'legendonly',
            line:{color:'#8ecae6',width:1.2}
          },
		  {
		    x: <?=json_encode($bpd)?>,
		    y: <?=json_encode($bpv)?>,
		    name:'BTC 12M Change %',
		    yaxis:'y2',
		    line:{color:'#FF0000'},
		    customdata: <?=json_encode($btc_prices_hover_plot)?>,
		    hovertemplate:
		      '<b>BTC 12M Change</b>: %{y:.2f}%<br>' +
		      '<b>BTC Price</b>: $%{customdata:,.0f}<br>' +
		      '<extra></extra>'
		  },
          {
            x: <?=json_encode($sp500_12m_dates)?>,
            y: <?=json_encode($sp500_12m_vals)?>,
            name:'SP500 12M Change %',
            yaxis:'y5',
            visible:'legendonly',
            line:{color:'#ff2bd6',width:2.2},
            customdata: <?=json_encode($sp500_12m_levels)?>,
            hovertemplate:
              '<b>SP500 12M Change</b>: %{y:.2f}%<br>' +
              '<b>SP500</b>: %{customdata:,.2f}<br>' +
              '<extra></extra>'
          },
          {
            x: <?=json_encode($nasdaq_12m_dates)?>,
            y: <?=json_encode($nasdaq_12m_vals)?>,
            name:'NASDAQ 12M Change %',
            yaxis:'y6',
            visible:'legendonly',
            line:{color:'#b197fc',width:1.6},
            customdata: <?=json_encode($nasdaq_12m_levels)?>,
            hovertemplate:
              '<b>NASDAQ 12M Change</b>: %{y:.2f}%<br>' +
              '<b>NASDAQ</b>: %{customdata:,.2f}<br>' +
              '<extra></extra>'
          },
          {
            x: <?=json_encode($russell_12m_dates)?>,
            y: <?=json_encode($russell_12m_vals)?>,
            name:'Russell 2000 12M Change % (+5M)',
            yaxis:'y7',
            showlegend:true,
            visible:'legendonly',
            line:{color:'#ff922b',width:1.6},
            customdata: <?=json_encode($russell_12m_levels)?>,
            hovertemplate:
              '<b>Russell 2000 12M Change (+5M)</b>: %{y:.2f}%<br>' +
              '<b>Russell 2000</b>: %{customdata:,.2f}<br>' +
              '<extra></extra>'
          },
		  {
		    x: <?=json_encode($btc_peak_marker_dates)?>,
		    y: <?=json_encode($btc_peak_marker_vals)?>,
		    type:'scatter',
		    mode:'markers+text',
		    text: <?=json_encode(array_fill(0, count($btc_peak_marker_dates), 'Peak'))?>,
		    textposition:'top center',
		    yaxis:'y2',
		    name:'BTC Cycle Peaks',
		    marker:{
		      size:12,
		      symbol:'diamond',
		      color:'#ffd166',
		      line:{color:'#ffffff', width:1}
		    },
		    hovertemplate:'<b>BTC Cycle Peak</b><br>%{x}<br>BTC 12M Change: %{y:.2f}%<extra></extra>'
		  },
		  {
		    x: <?=json_encode($btc_bottom_marker_dates)?>,
		    y: <?=json_encode($btc_bottom_marker_vals)?>,
		    type:'scatter',
		    mode:'markers+text',
		    text: <?=json_encode(array_fill(0, count($btc_bottom_marker_dates), 'Bottom'))?>,
		    textposition:'bottom center',
		    yaxis:'y2',
		    name:'BTC Cycle Bottoms',
		    marker:{
		      size:12,
		      symbol:'diamond',
		      color:'#00d084',
		      line:{color:'#ffffff', width:1}
		    },
		    hovertemplate:'<b>BTC Cycle Bottom</b><br>%{x}<br>BTC 12M Change: %{y:.2f}%<extra></extra>'
		  }
		];
		
		Plotly.newPlot('chart', traces, {
		  paper_bgcolor:'#212529',
		  plot_bgcolor:'#212529',
		  font:{ color:'#fff' },
		  shapes: dominanceShapes,
		  annotations: [
		    {
		      xref:'paper',
		      yref:'paper',
		      x:0.01,
		      y:0.995,
		      xanchor:'left',
		      yanchor:'top',
		      align:'left',
		      showarrow:false,
		      bgcolor:'rgba(0,0,0,0.45)',
		      bordercolor:'rgba(255,255,255,0.12)',
		      borderwidth:1,
		      font:{size:12,color:'#fff'},
		      text:
		        '<b>Impact winner:</b> <?=htmlspecialchars($impact_winner)?>' +
		        '<br><b>Peak driver:</b> <?=htmlspecialchars($peak_driver_winner)?>' +
		        '<br><b>Bottom driver:</b> <?=htmlspecialchars($bottom_driver_winner)?>' +
		        '<br><b>Risk-Love impact:</b> <?=($risklove_impact_score!==null?number_format($risklove_impact_score,2):'n/a')?>' +
		        '<br><b>Current peak leader:</b> <?=htmlspecialchars($current_peak_condition_winner)?>' +
		        '<br><b>Current bottom leader:</b> <?=htmlspecialchars($current_bottom_condition_winner)?>' +
                '<br><b>Debt-liq latest gap:</b> <?=($debt_liq_latest_gap!==null?'$'.number_format($debt_liq_latest_gap,2).'T':'n/a')?>'
		    }
		  ],
		  hovermode:'x unified',
		  xaxis:{showgrid:false},
		  yaxis:{domain:[0.28,1], title:'Macro / Risk Waves', zeroline:false},
		  yaxis2:{overlaying:'y', side:'right', title:'BTC 12M Change %', showgrid:false, zeroline:false},
		  yaxis3:{domain:[0,0.20], title:'MOVE Index', showgrid:false, zeroline:false},
          yaxis4:{overlaying:'y', side:'left', position:0.04, title:'Debt Liquidity $T', showgrid:false, zeroline:false},
          yaxis5:{overlaying:'y', side:'right', position:0.94, title:'SP500 12M Change %', showgrid:false, zeroline:false},
          yaxis6:{overlaying:'y', side:'right', position:0.98, title:'NASDAQ 12M Change %', showgrid:false, zeroline:false},
          yaxis7:{overlaying:'y', side:'left', position:0.08, title:'Russell 2000 12M Change % (+5M)', showgrid:false, zeroline:false},
		  legend:{
		    orientation:'h',
		    y:-0.2,
		    x:0,
		    groupclick:'togglegroup',
		    itemclick:'toggle',
		    itemdoubleclick:'toggleothers'
		  }
		}).then(() => {
		  const gd = document.getElementById('chart');
		  const baseX = gd.data.map(t => (t.x ? t.x.slice() : []));
		  const shiftCache = new Map();
		
		  function shiftDatesCached(baseDates, months) {
		    const key = JSON.stringify([baseDates.length ? baseDates[0] : '', baseDates.length ? baseDates[baseDates.length-1] : '', months, baseDates.length]);
		    if (shiftCache.has(key)) return shiftCache.get(key);
		
		    const shifted = baseDates.map(d => {
		      const x = new Date(d + 'T00:00:00Z');
		      const day = x.getUTCDate();
		      x.setUTCMonth(x.getUTCMonth() + months);
		      if (x.getUTCDate() !== day) x.setUTCDate(0);
		      return x.toISOString().slice(0, 10);
		    });
		
		    shiftCache.set(key, shifted);
		    return shifted;
		  }
		
		  const BTC_SHIFT_MONTHS = 2;
		  const btcTraceIndex = gd.data.findIndex(t => t.name === 'BTC 12M Change %');
		  const peakTraceIndex = gd.data.findIndex(t => t.name === 'BTC Cycle Peaks');
		  const bottomTraceIndex = gd.data.findIndex(t => t.name === 'BTC Cycle Bottoms');
		
		  if (btcTraceIndex !== -1) {
		    const shiftedBTC = shiftDatesCached(baseX[btcTraceIndex], BTC_SHIFT_MONTHS);
		    Plotly.restyle(gd, { x: [shiftedBTC] }, [btcTraceIndex]);
		  }
		  if (peakTraceIndex !== -1) {
		    const shiftedPeaks = shiftDatesCached(baseX[peakTraceIndex], 0);
		    Plotly.restyle(gd, { x: [shiftedPeaks] }, [peakTraceIndex]);
		  }
		  if (bottomTraceIndex !== -1) {
		    const shiftedBottoms = shiftDatesCached(baseX[bottomTraceIndex], 0);
		    Plotly.restyle(gd, { x: [shiftedBottoms] }, [bottomTraceIndex]);
		  }
		
		  const offsets = {
		    macro: 37,
		    liq: -36,
		    ycc: -60,
		    qe: -27,
		    blend: -27,
		    risklove: -27,
            debtliq: -28
		  };
		
		  const updatesX = [];
		  const traceIndices = [];
		
		  gd.data.forEach((t, i) => {
		    if (!t.legendgroup) return;
		    const g = t.legendgroup;
		    if (!(g in offsets)) return;
		
		    const totalShift = offsets.macro + offsets[g];
		    updatesX.push(shiftDatesCached(baseX[i], totalShift));
		    traceIndices.push(i);
		  });
		
		  if (traceIndices.length) {
		    Plotly.restyle(gd, { x: updatesX }, traceIndices);
		  }
		
		  
          // Force Stealth QE deviation fills to render on initial load.
          // Without this post-shift restyle, Plotly can leave the tonexty fills unpainted
          // until the QE legend group is toggled off/on.
          const qeDeviationTraceIndices = [];
          gd.data.forEach((t, i) => {
            if (t.legendgroup === 'qe' && t.showlegend === false) qeDeviationTraceIndices.push(i);
          });
          if (qeDeviationTraceIndices.length) {
            Plotly.restyle(gd, { visible: true }, qeDeviationTraceIndices);
          }
        
          const currentEndDate = gd._fullLayout.xaxis.range[1];
		  Plotly.relayout(gd, {
		    'xaxis.range': ['2019-08-01', currentEndDate]
		  });
		}).catch(err => {
		  console.error('Plotly render error:', err);
		});
	</script>
	
	<script type="text/javascript">
		var contentHeight = $('body').outerHeight(true);
		if (window.frameElement) {
		  $(window.frameElement).height(contentHeight);
		}
	</script>
</body>
</html>
