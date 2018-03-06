import React from 'react';
import {Table, Column, Cell} from 'fixed-data-table-2';
import "fixed-data-table-2/dist/fixed-data-table.css";


class MyTextCell extends React.Component {

    render() {
        const {rowIndex, data, field} = this.props;
        return (
            <Cell {...this.props}>
                {data[rowIndex][field]}
            </Cell>
       );
    }
}


class MyTable extends React.Component {

    constructor() {
        super();

        this.state = {
            rows: [],
            filteredDataList: [],
            sortBy: 'Sample',
            sortDir: 'ASC',
            wheight: '0',
            wwidth: '0'
        };
        this.updateWindowDimensions = this.updateWindowDimensions.bind(this);
    }

    componentWillMount() {
        this.jsonList();
    }
    componentDidMount() {
        this.updateWindowDimensions();
        window.addEventListener('resize', this.updateWindowDimensions);
    }
    componentWillUnmount() {
        window.removeEventListener('resize', this.updateWindowDimensions);
    }
    updateWindowDimensions() {
        this.setState({ wwidth: window.innerWidth, wheight: window.innerHeight });
    }

    jsonList() {
        fetch('/content/metadata/').then(
            function(response){
                return response.json();
            }
        ).then(data => {
            var d = data;
            this.setState({
                rows: d.data,
                filteredDataList: d.data
            })
        })
    }

    render() {
        var sortDirArrow = '';
        if (this.state.sortDir !== null){
            sortDirArrow = this.state.sortDir === 'DESC' ? ' ↓' : ' ↑';
        }
        return <Table
            height={this.state.wheight - 150}
            width={this.state.wwidth - 100}
            rowsCount={this.state.filteredDataList.length}
            rowHeight={30}
            headerHeight={90}
            rowGetter={function(rowIndex) {return this.state.filteredDataList[rowIndex]; }.bind(this)}>
            <Column dataKey="cell_name" width={250} label={'Cell Name'+ (this.state.sortBy === 'cell_name' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="dataset" width={175} label={'Dataset' + (this.state.sortBy === 'dataset' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="cell_type" width={75} label={'Cell Type' + (this.state.sortBy === 'cell_type' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="global_mCH" width={100} label={'Global mCH' + (this.state.sortBy === 'global_mCH' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="global_mCG" width={100} label={'Global mCG' + (this.state.sortBy === 'global_mCG' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="global_mCA" width={100} label={'Global mCA' + (this.state.sortBy === 'global_mCA' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="global_mCCC" width={100} label={'Global mCCC' + (this.state.sortBy === 'global_mCCC' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="estimated_mCH" width={100} label={'Estimated mCH' + (this.state.sortBy === 'estimated_mCH' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="estimated_mCG" width={100} label={'Estimated mCG' + (this.state.sortBy === 'estimated_mCG' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="percent_genome_covered" width={100} label={'% Genome Covered' + (this.state.sortBy === 'percent_genome_covered' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="total_reads" width={100} label={'Total Reads' + (this.state.sortBy === 'total_reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="mapped_reads" width={100} label={'Mapped Reads' + (this.state.sortBy === 'mapped_reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="mapping_rate" width={100} label={'Mapping Rate' + (this.state.sortBy === 'mapping_rate' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="nonclonal_reads" width={100} label={'Nonclonal Reads' + (this.state.sortBy === 'nonclonal_reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="percent_nonclonal_rate" width={125} label={'% Nonclonal Rate' + (this.state.sortBy === 'percent_nonclonal_rate' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="filtered_reads" width={100} label={'Filtered Reads' + (this.state.sortBy === 'filtered_reads' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="filtered_rate" width={100} label={'Filtered Rate' + (this.state.sortBy === 'filtered_rate' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
            <Column  dataKey="lambda_mC" width={75} label={'Lambda mC' + (this.state.sortBy === 'lambda_mC' ? sortDirArrow : '')} headerRenderer={this._renderHeader.bind(this)}/>
        </Table>;
    }

    _onFilterChange(cellDataKey, event) {
        if (!event.target.value) {
            this.setState({
                filteredDataList: this.rows,
            });
        }
        var filterBy = event.target.value.toString().toLowerCase();
        var size = this.rows.length;
        var filteredList = [];
        for (var index = 0; index < size; index++) {
            var v = this.rows[index][cellDataKey];
            if (v.toString().toLowerCase().indexOf(filterBy) !== -1) {
                filteredList.push(this.rows[index]);
            }
        }
        this.setState({
            filteredDataList: filteredList,
        });
    }
    _sortRowsBy(cellDataKey) {
        var sortDir = this.state.sortDir;
        var sortBy = cellDataKey;
        if (sortBy === this.state.sortBy) {
            sortDir = this.state.sortDir === 'ASC' ? 'DESC' : 'ASC';
        } else {
            sortDir = 'DESC';
        }
        var rows = this.state.filteredDataList.slice();
        rows.sort((a, b) => {
            var sortVal = 0;
            if (a[sortBy] > b[sortBy]) {
                sortVal = 1;
            }
            if (a[sortBy] < b[sortBy]) {
                sortVal = -1;
            }

            if (sortDir === 'DESC') {
                sortVal = sortVal * -1;
            }
            return sortVal;
        });

        this.setState({sortBy, sortDir, filteredDataList : rows});
    }

    _renderHeader(label, cellDataKey) {
        return <div>
            <a onClick={this._sortRowsBy.bind(this, cellDataKey)}>{label}</a>
            <div>
                <input type="text" style={{width:90+'%'}} onChange={this._onFilterChange.bind(this, cellDataKey)}/>
            </div>
        </div>;
    }
}

module.exports = MyTable;
